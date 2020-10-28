# Scattering
## Installing the Scattering Package

Due to how Julia packages and Git interact, the Julia code is separated out into another package.
The package is available [here](https://github.com/ErlendLima/Scattering). For completeness' sake, its
installation instructions are repeated here.

Clone the package with
```shell
git clone https://github.com/ErlendLima/Scattering.git
```
or equivalent command for SSH.

Once cloned, the `Scattering` directory needs to be added to Julia's path. 
For Linux, the directory needs to be linked to `$HOME/.julia/dev`:
```shell
ln -s <PATH TO SCATTERING> $HOME/.julia/dev/
```
e.g.
```shell
ln -s ~/gits/Scattering/ $HOME/.julia/dev/
```
Finally the package must be added. Go to the `$HOME/.julia/dev/`, open Julia with the `julia` command, 
type `]` to enter package mode, and add the package with `add .Scattering`. E.g.
```shell
> cd ~/.julia/dev/
> julia
(@v1.5) pkg> add Scattering
```
The package can now be imported as any other package as `using Scattering`. Both
notebooks will work correctly and reproduce all results.
