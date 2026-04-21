# SFTP via Terminal

The `sftp` command-line client is available on Linux, macOS, and Windows
(via PowerShell or MSYS2).

## Prerequisites

- An SSH key registered with HIE-Storage. See [SSH Key Management](ssh-keys.md)
  and [HIE-Storage](hie-storage.md) for details.

## Connecting

```bash
sftp u123456@hie-storage.sstars.ws
```

Replace `u123456` with your actual username. If your private key is not at the
default location (`~/.ssh/id_ed25519`), specify it explicitly:

```bash
sftp -i /path/to/private_key u123456@hie-storage.sstars.ws
```

## Common Commands

Once connected, you are in an interactive session with two working directories:
the **remote** (on HIE-Storage) and the **local** (on your machine).

| Command | Description |
|---|---|
| `ls` | List remote directory contents |
| `lls` | List local directory contents |
| `cd <path>` | Change remote directory |
| `lcd <path>` | Change local directory |
| `pwd` | Print remote working directory |
| `lpwd` | Print local working directory |
| `get <file>` | Download a single file |
| `get -r <dir>` | Download a directory recursively |
| `mget <pattern>` | Download files matching a glob pattern |
| `put <file>` | Upload a single file |
| `exit` | Close the connection |

## Example: Downloading OzFlux Forcing Data

```bash
sftp u123456@hie-storage.sstars.ws
sftp> cd /data/hiestorage/WorkingData/MEDLYN_GROUP/PROJECTS/dynamics_simulations/lpj-guess/forcing/ozflux
sftp> lcd /path/to/lpj-guess-dave/benchmarks/ozflux
sftp> get -r processed forcing
sftp> exit
```

This downloads the `processed` directory and saves it locally as `forcing`.

## Non-Interactive Download with scp

For scripted or single-command transfers, `scp -r` is often more convenient
than `sftp`:

```bash
scp -r u123456@hie-storage.sstars.ws:/path/to/remote/dir /local/destination/
```
