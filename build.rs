// For embedding application icons, in Windows.

use std::{env, io};

use winresource::WindowsResource;

fn main() -> io::Result<()> {
    if env::var_os("CARGO_CFG_WINDOWS").is_some() {
        WindowsResource::new()
            // This path can be absolute, or relative to your crate root.
            .set_icon("src/resources/icon.ico")
            .compile()?;
    }
    Ok(())
}
