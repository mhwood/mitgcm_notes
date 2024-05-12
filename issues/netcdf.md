I ran into an issue with netcdf so I am adding this note

If your compilation can't find the netcdf path, check with this:
```
nc-config --all
```
It will print a lot of info including the paths to your include and lib:
```
  --prefix        -> /opt/homebrew/Cellar/netcdf/4.9.2_1
  --includedir    -> /opt/homebrew/Cellar/netcdf/4.9.2_1/include
  --libdir        -> /opt/homebrew/Cellar/netcdf/4.9.2_1/lib
```
