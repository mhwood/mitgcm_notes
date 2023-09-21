## Changing Permissions

This isn't really an MITgcm note, but it pertains to permissions and sharing files on a Unix system.

If you would like to give full permissions to another user on the systenm try
```
setfacl -R -m u:[user_name]:rwx [path_to_directory]
```

