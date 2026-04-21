# SFTP via FileZilla

[FileZilla](https://filezilla-project.org/) is a graphical SFTP client
available on Windows, macOS, and Linux.

## Prerequisites

- FileZilla installed (download from [filezilla-project.org](https://filezilla-project.org/))
- An SSH key registered with HIE-Storage. See [SSH Key Management](ssh-keys.md)
  and [HIE-Storage](hie-storage.md) for details.

## Connecting to HIE-Storage

Open the Site Manager (top-left button, or `File → Site Manager`). Click
**New site** and enter the following:

| Field | Value |
|---|---|
| Protocol | SFTP – SSH File Transfer Protocol |
| Host | `hie-storage.sstars.ws` |
| Port | *(leave empty)* |
| Logon Type | Key file |
| User | Your HIE-Storage username (e.g. `u123456`) |
| Key file | Path to your **private** key (e.g. `~/.ssh/id_ed25519`) |

Click **OK** to save. To connect, select the site from the dropdown to the
right of the Site Manager button and click **Connect**.

!!! note
    Always select the private key file (no `.pub` extension), not the public
    key. The public key (`.pub` file) is what you register with servers; the
    private key is what you use to authenticate.

## Downloading Files

Navigate to the remote directory on the right-hand panel and the local
destination on the left-hand panel. To download:

- **Single file or folder**: right-click → **Download**
- **Multiple items**: select them (Shift+click or Ctrl+click), then right-click
  → **Download**

Large recursive downloads can be queued and will transfer in the background.
