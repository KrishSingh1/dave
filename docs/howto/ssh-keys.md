# SSH Key Management

SSH keys are used for:

- Authenticating with [HIE-Storage](hie-storage.md) (SFTP)
- Authenticating with the [Lund GitLab](https://stormbringer4.nateko.lu.se/)
  server for `git clone`/push over SSH

The same key pair can be used for both.

## Key Types

**ed25519** is the recommended key type — it is smaller, faster, and more
secure than RSA. HIE-Storage explicitly rejects RSA keys.

## Checking for an Existing Key

Before generating a new key, check whether you already have one:

**Linux / macOS**
```bash
ls ~/.ssh/
```

**Windows**

Navigate to `C:\Users\<username>\.ssh\` in Explorer, or run:
```
dir C:\Users\<username>\.ssh\
```

Look for a file pair such as `id_ed25519` (private key) and `id_ed25519.pub`
(public key). If these exist, do not overwrite them unless you know what you
are doing — you would need to re-register the new key everywhere the old one
is used.

## Generating a New Key

```bash
ssh-keygen -t ed25519
```

Accept the default file location by pressing Enter. You may optionally set a
passphrase; if you leave it blank, the key will not require a password on use.

The command creates two files:

| File | Description |
|---|---|
| `~/.ssh/id_ed25519` | **Private key** — treat like a password, never share |
| `~/.ssh/id_ed25519.pub` | **Public key** — safe to share; registered with services |

## Registering the Key with GitLab

1. Copy the contents of your public key (id_ed25519.pub):

This can be opened in a text editor, or printed from a terminal:

**Windows**

```cmd
type %USERPROFILE%\.ssh\id_ed25519.pub
```

**Linux/MacOS**

```bash
cat ~/.ssh/id_ed25519.pub
```

2. In GitLab, go to your
   [preferences page](https://stormbringer4.nateko.lu.se/-/profile/preferences),
   click **SSH Keys** in the left sidebar, then **Add new key**.

3. Paste the public key content and click **Add key**.

## Registering the Key with HIE-Storage

Contact Drew to have your public key (`id_ed25519.pub`) added to your
HIE-Storage account. Do not send the private key.
