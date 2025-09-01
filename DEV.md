# Developer notes

## Artifacts

### What I did

1. Download `76Ge_0vbb.zip` from `https://nucleartheory.yale.edu/sites/default/files/files/76Ge_0vbb.zip`.
2. Run `unzip 76Ge_0vbb.zip`  to get relevant files.
3. Run `mkdir 76Ge_0vbb && mv *.txt 76Ge_0vbb` to organize the files in the directory.
1. Run `git init && git add 76Ge_0vbb && git write-tree` to get the tree hash of the directory (`762f23229b4be70e6e5a65744cd63e6e83eae8b4`).
1. Run `rm -rf .git` to clean up the git repository, as we only need the tree hash.
1. Run `tar -czf 76Ge_0vbb.tar.gz 76Ge_0vbb` to create a compressed tarball of the directory.
1. Run `sha256sum 76Ge_0vbb.tar.gz` to get the SHA256 checksum of the tarball (`83bf36e9a0ed23c8c97a67cd8b9a7fd76cb0f1ccdf328afca25aa77831e411be`).
1. Arrange the above information in `Artifacts.toml`.
3. Repeat for `76Ge_0vbb`

### Remarks

I tried to do the above on my MacBook Pro first but macOS added some hidden files like `._76Ge_0vbb` to the decompressed directory, which invalidates the tree hash. I moved to NERSC and everything works fine.
