//! The About popup: what this program is, which version is running, whether there's a newer one,
//! where to find it online, and where it keeps its preferences.

use std::{env, io, path::Path, process::Command, time::Duration};

use eframe::egui::{
    Color32, Context, CursorIcon, Label, ProgressBar, RichText, Sense, Ui, ViewportCommand, Window,
};
use software_manager::{Application, UpdateStatus};

use crate::{
    VERSION,
    gui::{ROW_SPACING, theme::COLOR_INFO},
    state::State,
};

const REPO_URL: &str = "https://github.com/David-OConnor/plascad";
const HOME_URL: &str = "https://www.athanorlab.com/plascad";

/// A text hyperlink: underlined, and opens the user's browser when clicked.
fn link(ui: &mut Ui, text: &str, url: &str) {
    let resp = ui
        .add(Label::new(RichText::new(text).color(COLOR_INFO).underline()).sense(Sense::click()))
        .on_hover_cursor(CursorIcon::PointingHand);

    if resp.clicked()
        && let Err(e) = webbrowser::open(url)
    {
        eprintln!("Failed to open the web browser: {:?}", e);
    }
}

/// Open a folder in the OS's file browser.
fn open_dir(dir: &Path) -> io::Result<()> {
    let program = if cfg!(windows) {
        "explorer"
    } else if cfg!(target_os = "macos") {
        "open"
    } else {
        "xdg-open"
    };

    Command::new(program).arg(dir).spawn().map(|_| ())
}

/// Check GitHub for a newer release, download it, and install it. Installing waits for PlasCAD to
/// close, since its files can't be replaced while it runs; the release's own setup script then
/// installs it, and PlasCAD reopens.
fn update_section(state: &mut State, ui: &mut Ui) {
    let updater = &mut state.ui.updater;
    updater.poll();

    if updater.busy() {
        // So the download progress moves.
        ui.ctx().request_repaint_after(Duration::from_millis(100));
    }

    let mut check = false;
    let mut download = false;
    let mut install = false;

    ui.horizontal(|ui| match updater.status() {
        UpdateStatus::Idle => {
            check = ui.button("Check for updates").clicked();
        }
        UpdateStatus::Checking(_) => {
            ui.spinner();
            ui.label("Checking for updates…");
        }
        UpdateStatus::UpToDate => {
            ui.label(RichText::new("This is the latest version.").color(Color32::WHITE));
            check = ui.button("Check again").clicked();
        }
        UpdateStatus::Available(release) => {
            ui.label(
                RichText::new(format!("Version {} is available.", release.version))
                    .color(Color32::WHITE),
            );
            link(ui, "Release notes", &release.page_url);

            download = ui
                .button("Download and install")
                .on_hover_text(format!(
                    "Download the update ({:.0} MB). Installing it is the next step.",
                    release.size as f64 / 1_000_000.
                ))
                .clicked();
        }
        UpdateStatus::Downloading {
            release, progress, ..
        } => {
            ui.add(
                ProgressBar::new(progress.fraction())
                    .text(format!("Downloading {}…", release.version))
                    .desired_width(300.),
            );
        }
        UpdateStatus::Ready(staged) => {
            ui.label(format!("Version {} is downloaded.", staged.release.version));

            install = ui
                .button("Restart and install")
                .on_hover_text(
                    "Close PlasCAD, install the update, then reopen it. Preferences are saved \
                    first; save your files before this.",
                )
                .clicked();
        }
        UpdateStatus::Installing => {
            ui.label("Installing…");
        }
        UpdateStatus::Failed(e) => {
            ui.label(
                RichText::new(format!("Unable to update: {e}"))
                    .color(Color32::ORANGE)
                    .small(),
            );
            check = ui.button("Try again").clicked();
        }
    });

    let ctx = ui.ctx().clone();
    let repaint = {
        let ctx = ctx.clone();
        move || ctx.request_repaint()
    };

    if check {
        updater.check(Application::PlasCad, VERSION, repaint);
    } else if download {
        updater.download(repaint);
    } else if install && updater.install_on_exit() {
        // The installer is waiting for this process to end.
        state.update_save_prefs();
        ctx.send_viewport_cmd(ViewportCommand::Close);
    }
}

fn contents(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.heading(RichText::new("PlasCAD").color(Color32::WHITE));
        ui.label(RichText::new(format!("v{VERSION}")).color(Color32::GRAY));
    });

    ui.add_space(ROW_SPACING / 2.);

    ui.label("Plasmid and primer design, PCR, and related tools.");

    ui.add_space(ROW_SPACING / 2.);

    update_section(state, ui);

    ui.add_space(ROW_SPACING / 2.);

    ui.horizontal(|ui| {
        ui.label("Home page:");
        link(ui, HOME_URL, HOME_URL);
    });

    ui.horizontal(|ui| {
        ui.label("Source code:");
        link(ui, REPO_URL, REPO_URL);
    });

    ui.add_space(ROW_SPACING / 2.);

    // The preferences file is relative to the working directory; the Start menu entry and desktop
    // file start PlasCAD in its install folder.
    ui.label("Preferences are saved in the folder PlasCAD was started from:");

    let dir = env::current_dir().ok();

    ui.horizontal(|ui| {
        if let Some(dir) = &dir {
            if ui.button("Open folder").clicked() {
                state.ui.about_folder_error = open_dir(dir).err().map(|e| e.to_string());
            }

            ui.label(RichText::new(dir.display().to_string()).color(Color32::WHITE));
        }
    });

    if let Some(error) = &state.ui.about_folder_error {
        ui.label(
            RichText::new(format!("Could not open the folder: {error}"))
                .color(Color32::ORANGE)
                .small(),
        );
    }
}

pub fn about_window(state: &mut State, ctx: &Context) {
    if !state.ui.about {
        return;
    }

    // Kept separate from the state while the contents borrow it.
    let mut open = true;

    Window::new("About")
        .open(&mut open)
        .default_width(480.)
        .resizable(true)
        .collapsible(false)
        .show(ctx, |ui| contents(state, ui));

    state.ui.about &= open;
}
