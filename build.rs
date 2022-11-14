extern crate bindgen;

// use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let _output = Command::new("make")
        .current_dir("c-impl")
        .arg("ARCH=1")
        .output()
        .expect("make command failed");

    // let profile = std::env::var("PROFILE").expect("Failed to get build profile");
    let out_path = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    // let target = std::env::var("TARGET").expect("Failed to get target");

    // Tell cargo to look for shared libraries in the specified directory
    println!("cargo:rustc-link-search=c-impl/PQCrypto-SIDH/lib434poik",);

    // Tell cargo to tell rustc to link the sidh library.
    println!("cargo:rustc-link-lib=sidh");

    // Tell cargo to invalidate the built crate whenever the wrapper changes
    println!("cargo:rerun-if-changed=c-impl/bindgen-wrapper.h");

    let bindings = bindgen::Builder::default()
        .header("c-impl/bindgen-wrapper.h")
        .derive_debug(true)
        .derive_default(true)
        .derive_eq(true)
        .derive_hash(true)
        .derive_partialeq(true)
        .generate()
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
