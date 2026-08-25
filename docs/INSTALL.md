# Installing NeoEdit

**Download page: https://github.com/evozoa/NeoEdit/releases/latest**

Pick the file for your computer, then follow the one-time steps below. No other software is needed.

| Computer | Download |
|---|---|
| **Windows** (10 / 11) | [NeoEdit-Windows-Setup.exe](https://github.com/evozoa/NeoEdit/releases/latest/download/NeoEdit-Windows-Setup.exe) |
| **Mac with Apple Silicon** (M1, M2, M3, M4 — most Macs from 2021 on) | [NeoEdit-macOS-AppleSilicon.dmg](https://github.com/evozoa/NeoEdit/releases/latest/download/NeoEdit-macOS-AppleSilicon.dmg) |
| **Mac with an Intel chip** (older Macs) | [NeoEdit-macOS-Intel.dmg](https://github.com/evozoa/NeoEdit/releases/latest/download/NeoEdit-macOS-Intel.dmg) |

Not sure which Mac you have? Click the Apple menu () ▸ **About This Mac**. The *Chip* line says
*Apple M1/M2/M3/M4* (→ Apple Silicon) or *Intel* (→ Intel).

## Windows

1. Open the downloaded **NeoEdit-Windows-Setup.exe** (usually in your *Downloads* folder).
2. If a blue box says **"Windows protected your PC"**, click **More info**, then **Run anyway**.
   NeoEdit is open-source but not code-signed by a company, so Windows shows this warning once.
3. Click **Next** through the installer. It installs for your user account, so no administrator
   password is needed. You can tick "Create a desktop icon".
4. Find **NeoEdit** in the Start Menu (or on the desktop).

To remove it later: Settings ▸ Apps ▸ NeoEdit ▸ Uninstall.

## Mac

1. Open the downloaded **.dmg** file. A window appears with NeoEdit and an *Applications* folder:
   drag **NeoEdit** onto **Applications**.
2. Open your **Applications** folder and double-click **NeoEdit**. macOS says
   *"Apple could not verify NeoEdit is free of malware…"* — click **Done** (not "Move to Trash").
3. Open **System Settings ▸ Privacy & Security**, scroll down to the **Security** section, and
   click **Open Anyway** next to the NeoEdit message. Enter your Mac password and click **Open**.
4. Done — from now on NeoEdit opens normally from Applications, Launchpad or Spotlight.

You only do steps 2–3 once. On macOS 13 or 14 you can instead right-click NeoEdit ▸ **Open** ▸ **Open**.

To remove it later: drag NeoEdit from Applications to the Trash.

## What is included

* NeoEdit itself, with the example data under *File ▸ Open recent* once you have opened something.
* **MAFFT** for *Alignment ▸ MAFFT* — no separate install.
* Everything works offline except **File ▸ Import from NCBI / Ensembl / UCSC**, which needs
  an internet connection.

## Trouble?

* **Windows: "This app can't run on your PC"** — you downloaded the Mac file; use the `.exe`.
* **Mac: "NeoEdit is damaged and can't be opened"** — the download was interrupted or the
  wrong chip version was chosen. Delete it, download the other `.dmg`, and repeat the steps.
* **Mac: the "Open Anyway" button is missing** — double-click NeoEdit first (step 2); the
  button appears in Privacy & Security only after that attempt, and only for a few minutes.
* Anything else: https://github.com/evozoa/NeoEdit/issues
