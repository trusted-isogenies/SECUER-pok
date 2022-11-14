//! Generate, prove, and verify isogeny walks as part of a SECUER trusted setup.
//!
//! ```rust
//! let secret = Secret::generate(*J_0).expect("secret generation usually should suceed");
//!
//! let mut proof =
//!     Proof::try_from(secret).expect("proving knowledge of the secret should succeed");
//!
//! assert!(proof.verify().is_ok())
//! ```

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

use std::{
    collections::HashSet,
    fmt::{self, Debug},
    hash::{Hash, Hasher},
};

use hex;
use lazy_static::lazy_static;
use subtle::ConstantTimeEq;
use thiserror::Error;

// mod serialize;

pub mod bindings {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

lazy_static! {
    /// j-invariant = 0
    ///
    /// Identifier (up to isomorphism) of the short Weierstrass elliptic curve
    /// E₀: y² = x³ + x, as defined over 𝔽p². A 'special' (for some definition
    /// of special) curve that is the common starting curve for isogeny
    /// protocols.
    ///
    /// For the generation of a [SECUER] curve, this will be the first curve in
    /// the isogeny chain, and never repeated.
    ///
    /// [SECUER]: https://eprint.iacr.org/2022/1469.pdf
    static ref J_0: JInvariant = Default::default();
}

/// A j-invariant defined over 𝔽p² that uniquely defines an isomorphism class of
/// elliptic curves.
///
/// The j-invariant functions as an identifier of an elliptic curve, up to
/// isomorphism.
#[derive(Clone, Copy, Default)]
pub struct JInvariant(bindings::packed_fp2);

impl Debug for JInvariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // f.debug_struct("JInvariant")
        // .field("real", &self.0.real)
        // .field("imaginary", &self.0.imag)
        // .finish()
        unsafe {
            f.debug_tuple("JInvariant")
                .field(&hex::encode(&self.0.bytes))
                .finish()
        }
    }
}

impl Eq for JInvariant {}

impl Hash for JInvariant {
    fn hash<H: Hasher>(&self, state: &mut H) {
        unsafe { state.write(&self.0.bytes) }
    }
}

impl PartialEq for JInvariant {
    fn eq(&self, other: &Self) -> bool {
        unsafe { self.0.bytes.ct_eq(&other.0.bytes).into() }
    }
}

// impl serialize::Serialize for JInvariant {}

// impl serialize::Deserialize JInvariant {}

/// A secret isogeny walk from a domain curve `E` uniquely specified by the
/// kernel generator `P`.
///
/// `E` is comprised of `A` and `C`, elements of 𝔽p², used to compute the
/// canconical j-invariant of the Montgomery representation of the curve.
///
/// `P` is a point defined on `E` over 𝔽p² that generates the finite subgroup of
/// points comprising the kernel of the secret isogeny.
///
/// This walk in the supersingular isogeny graph is what we [prove] knowledge
/// of.
///
/// [prove]: https://eprint.iacr.org/2022/1469.pdf
#[derive(Default)]
pub struct Secret(bindings::secret);

impl Secret {
    /// Randomly generates a secret isogeny walk from the domain curve argument
    /// to a new supersingular curve.
    pub fn generate(domain_curve: &JInvariant) -> Result<Secret, Error> {
        let mut sec: bindings::secret = Default::default();

        unsafe {
            if bindings::Generate(&domain_curve.0, &mut sec) {
                Ok(Self(sec))
            } else {
                Err(Error::SecretGenerationFailure)
            }
        }
    }
}

/// A [proof of knowledge] that a cyclic separable isogeny of smooth degree
/// exists from a supersingular curve `E0` to `E1` over 𝔽p².
///
/// The `response` is comprised of a series of challenge trits (-1, 0, 1 for the
/// left, bottom, or right edges of an SIDH square respectively), and the
/// corresponding isogeny response (in the Σ-protocol sense).
///
/// [proof of knowledge]: https://eprint.iacr.org/2022/1469.pdf
#[derive(Clone, Default)]
pub struct Proof(Box<bindings::proof>);

impl Proof {
    pub fn verify(&mut self) -> Result<(), Error> {
        unsafe {
            if bindings::Verify(&mut *self.0) {
                Ok(())
            } else {
                Err(Error::VerificationFailure)
            }
        }
    }

    /// Get the domain (start) supersingular elliptic curve `E₀` as its j-invariant.
    pub fn domain(&self) -> JInvariant {
        JInvariant(self.0.E0)
    }

    /// Get the codomain (destination) supersingular elliptic curve `E₀` as its j-invariant.
    pub fn codomain(&self) -> JInvariant {
        JInvariant(self.0.E1)
    }
}

impl TryFrom<Secret> for Proof {
    type Error = Error;

    fn try_from(secret: Secret) -> Result<Proof, Self::Error> {
        let mut proof: Box<bindings::proof> = Default::default();

        unsafe {
            // TODO: this may panic; if it does, catch and return Error::ProvingFailure
            bindings::Prove(&secret.0, &mut *proof);
        }

        Ok(Proof(proof))
    }
}

// impl serialize::Serialize for Proof {}

// impl serialize::Deserialize for Proof {}

/// Proving errors
#[derive(Debug, Error)]
pub enum Error {
    #[error("secret isogeny chain generation failed")]
    SecretGenerationFailure,

    #[error("proving knowledge of the secret isogeny failed")]
    ProvingFailure,

    #[error("verifying the proof failed")]
    VerificationFailure,
}

/// A structure to track the series of isogeny walks conducted during SECUER trusted setup.
///
/// We want to keep the in-order hops from curve to curve and the proofs of knowledge of
/// the secret isogeny, and we want to keep a handy unique cached set of the observed curves.
///
/// There MUST be NO duplicated curves in the entire series, such a 'collision'
/// indicates a cycle in the complete secret isogeny walk.
#[derive(Clone)]
pub struct Sequence {
    /// All the curves that have been observed in this sequence of isogeny
    /// walks, including the start curve `E₀` and the final SECUER curve of
    /// unknown endomorphism ring.
    curves: HashSet<JInvariant>,
    /// The sequence of isogeny walks, tracked by the verified proofs of isogeny
    /// knowledge and their associated curves.
    hops: Vec<Proof>,
}

impl Sequence {
    /// Start a new [`Sequence`] at the `E₀` curve.
    pub fn new() -> Self {
        let mut sequence = Sequence {
            curves: std::collections::HashSet::new(),
            hops: Vec::new(),
        };

        sequence.curves.insert(*J_0);

        return sequence;
    }

    /// Get the current tip curve, the codomain of the last/most recent hop in the isogeny walk.
    pub fn tip_curve(&self) -> JInvariant {
        // If there is only one curve in our set and no hops, return that.
        if (self.curves.len() == 1) && (self.hops.len() == 0) {
            // There's only one curve this could be...
            return J_0.clone();
        } else {
            self.hops
                .last()
                .expect("there should be a tip with a codomain")
                .codomain()
                .clone()
        }
    }

    /// Add a new isogeny hop to the end of this SECUER sequence.
    ///
    /// Verifies the proof of knowledge of an isogeny between the tip curve in
    /// the current sequence and the codomain curve proposed to be added.
    ///
    /// Also verifies that the new codomain curve has never been observed in
    /// this sequence before.
    pub fn push(&mut self, proof: &mut Proof) -> Result<(), Error> {
        // TODO: stop panic'ing! Handle these errors appropriately

        // New isogeny hop MUST have our known tip curve as its domain curve for the isogeny walk.
        assert_eq!(self.tip_curve(), proof.domain());

        let new_tip_curve = proof.codomain();

        // The new tip curve MUST NOT have been seen before in this sequence.
        assert!(!self.curves.contains(&new_tip_curve));

        // Proof of knowledge of an isogeny between the tip curve and a new curve MUST verify
        assert!(proof.verify().is_ok());

        // OK! We can extend our sequence with this [`Proof`].
        self.curves.insert(new_tip_curve);
        self.hops.push(proof.clone());

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{
        bindings::{self, packed_fp2, proof, secret},
        *,
    };

    #[test]
    fn prove_and_verify_bindings() {
        unsafe {
            let startE: packed_fp2 = Default::default();
            let mut sec: secret = Default::default();

            assert!(bindings::Generate(&startE, &mut sec));
            // println!("secret: {:?}", sec.__bindgen_anon_1);

            let mut proof: proof = Default::default();

            bindings::Prove(&sec, &mut proof);
            // println!("proof: {:?}", stringify!(proof));

            assert!(bindings::Verify(&mut proof))
        }
    }

    #[test]
    fn prove_and_verify() {
        let secret = Secret::generate(&*J_0).expect("secret generation usually should succeed");

        // TODO: clone inside verify() to make mut?
        let mut proof =
            Proof::try_from(secret).expect("proving knowledge of the secret should succeed");

        assert!(proof.verify().is_ok())
    }

    #[test]
    fn single_party_secuer() {
        // let mut sequence: Vec<Proof> = Default::default();
        let mut sequence = Sequence::new();

        println!("Starting curve:");
        println!("{:?}", sequence.tip_curve());

        // Do a sequence of random isogeny walks in the supersingular curve
        // graph, as if each loop was an independent participant contibution.
        for i in 1..111 {
            println!("Embarking on hop {:?}", i);

            // Get the domain curve for a new isogeny walk from the codomain
            // curve of the current tip of the sequence.
            let domain_curve = sequence.tip_curve();

            let secret =
                Secret::generate(&domain_curve).expect("secret generation usually should succeed");

            // !!! stack overflow !!!
            let mut proof =
                Proof::try_from(secret).expect("proving knowledge of the secret should succeed");

            // Only submit if our proof succeeds.
            assert!(proof.verify().is_ok());

            // Hopefully it succeeds!
            assert!(sequence.push(&mut proof).is_ok());

            println!("Curve {:?}:", i);
            println!("{:?}", sequence.tip_curve());
        }

        println!("\n");
        println!("Final SECUER:");
        println!("{:?}", sequence.tip_curve())
    }
}
