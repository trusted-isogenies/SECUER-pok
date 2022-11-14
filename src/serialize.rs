//! Serialization and Deserialization traits
//!
//! We're basically defining our custom byte and hex serialization format, so we
//! define these traits ourselves instead of pulling in all of `serde`.

use std::io;

use thiserror::Error;

/// Serialization trait
pub trait Serialize: Sized {
    /// Write `self` to the given `writer` using the canonical format.
    ///
    /// Notice that the error type is [`std::io::Error`]; this indicates that
    /// serialization MUST be infallible up to errors in the underlying writer.
    /// In other words, any type implementing `Serialize` must make illegal
    /// states unrepresentable.
    fn serialize<W: io::Write>(&self, writer: W) -> Result<(), io::Error>;

    /// Helper function to construct a vec to serialize the current struct into
    fn serialize_to_vec(&self) -> Result<Vec<u8>, io::Error> {
        let mut data = Vec::new();
        self.serialize(&mut data)?;
        Ok(data)
    }
}

/// Deserialization trait
pub trait Deserialize: Sized {
    /// Try to read `self` from the given `reader`.
    fn deserialize<R: io::Read>(reader: R) -> Result<Self, SerializationError>;
}

/// A serialization error.
#[derive(Error, Debug)]
pub enum SerializationError {
    /// An io error that prevented deserialization
    #[error("io error: {0}")]
    Io(#[from] io::Error),

    /// The data to be deserialized was malformed.
    #[error("parse error: {0}")]
    Parse(&'static str),
}
