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

use lazy_static::lazy_static;
use thiserror::Error;

mod serialize;

mod bindings {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

/// A j-invariant defined over 𝔽p² that uniquely defines an isomorphism class of
/// elliptic curves.
///
/// A j-invariant functions as an identifier of an elliptic curve, up to
/// isomorphism.
#[derive(Clone, Copy, Default)]
pub struct JInvariant(bindings::packed_fp2);

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
    pub fn generate(domain_curve: JInvariant) -> Result<Secret, Error> {
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
#[derive(Clone, Copy, Default)]
pub struct Proof(bindings::proof);

impl Proof {
    pub fn verify(&mut self) -> Result<(), Error> {
        unsafe {
            if bindings::Verify(&mut self.0) {
                Ok(())
            } else {
                Err(Error::VerificationFailure)
            }
        }
    }
}

impl TryFrom<Secret> for Proof {
    type Error = Error;

    fn try_from(secret: Secret) -> Result<Proof, Self::Error> {
        let mut proof: bindings::proof = Default::default();

        unsafe {
            // TODO: this may panic; if it does, catch and return Error::ProvingFailure
            bindings::Prove(&secret.0, &mut proof);
        }

        Ok(Proof(proof))
    }
}

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
        let secret = Secret::generate(*J_0).expect("secret generation usually should suceed");

        // TODO: clone inside verify() to make mut?
        let mut proof =
            Proof::try_from(secret).expect("proving knowledge of the secret should succeed");

        assert!(proof.verify().is_ok())
    }
}
