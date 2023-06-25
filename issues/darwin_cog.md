# Darwin cog issues

The Darwin package uses a package called `cogapp` when compiling. I've run into 2 issues with this approach.

First, the `tools/darwin/cog` script contains a reference to `python`. Sometimes (e.g. on MacOS or an HPC if aliases/paths aren't set), then the line
```
#!/usr/bin/env python
```
will throw an error. To remedy, update `python` to `python3`.

Second, the `cogapp` class uses `hashlib.md5` which can throw the following error after recent HPC updates:
```
[digital envelope routines: EVP_DigestInit_ex] disabled for FIPS
```
To remedy this issue, I edited the tools/darwin/cogapp.py file so that calls to `md5` were passed the following argument:
```
hasher = hashlib.md5(usedforsecurity=False)
```
After this update, the model compiled as before.
