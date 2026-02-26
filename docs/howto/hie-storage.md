# HIE-Storage

HIE-Storage is a file server used by the DAVE team to store climate forcing
data and other large datasets.

## Connection Details

| Parameter | Value |
|---|---|
| Host | `hie-storage.sstars.ws` |
| Protocol | SFTP |
| Authentication | SSH key only (password authentication is disabled) |
| Username | `u<staff_id>` — e.g. `u123456` |

If you do not have an account, contact Drew.

## SSH Key Requirements

Password authentication is disabled. You must use an SSH key. **ed25519 keys
are recommended** — RSA keys are rejected as insecure.

If you do not have an SSH key yet, or need to create an ed25519 key, see
[SSH Key Management](ssh-keys.md).

## Accessing Files

There are two ways documented in this wiki to transfer files to and from
HIE-Storage:

- **[SFTP via FileZilla](sftp-filezilla.md)** — graphical client, suitable for
  browsing and downloading individual files or directories.
- **[SFTP via terminal](sftp-terminal.md)** — command-line access, useful for
  scripted transfers and working on remote systems.

## OzFlux Forcing Data

The climate forcing data required for OzFlux site simulations is at:

```
/data/hiestorage/WorkingData/MEDLYN_GROUP/PROJECTS/dynamics_simulations/lpj-guess/forcing/ozflux/processed
```

This directory should be downloaded and placed at
`benchmarks/ozflux/forcing` inside the repository. See
[OzFlux Site Simulations](ozflux-simulations.md) for the expected directory
structure.
