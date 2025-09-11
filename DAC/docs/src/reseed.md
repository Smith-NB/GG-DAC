# Reseeding

```@meta
CurrentModule = DAC
```
Reseeding is where the search restarts from a random starting structure periodically. There are currently two `Reseeder` types implemented. Any new types must flfull the following requirements.

An `updateHopsToReseed!` method must be implemented to track and update how many hops are left until a reseed is triggered.

A `checkNewlyAcceptedStructure!` method must be implemented so that, when a new structure is accepted, the `Reseeder` can check if this triggers any behaviour (e.g. if a new LES is found, this may reset how many hops are left until the next reseed occurs.)

A `timeToReseed!` function must be implemented. This function must return `true` or `false` depending on if it is time for a reseed to be performed.

A `getReseedPeriod` function must be implemented and should simply return the reseed period of the `Reseeder`.

The `Reseeder` must contain a field named `getReseedStructure` if `timeToReseed!` can return `true`. (If this sounds odd, note that the [`ReseedDisabled`](@ref) `Reseeder` will always return `false` on `timeToReseed!`). `getReseedStructure` is a `Function` passed to the `Reseeder` to generate a new random structure upon a reseed being triggered. Typically, `getReseedStructure` is set to [`generateRandomSeed`](@ref), an inbuilt function to return the structure of a random seed. Additionally, an `args` field must be defined for the `Reseeder` containing all the arguments needed for the `getReseedStructure` function. `args` will then be unpacked using the `...` operator for the `getReseedStructure` function. 

# Implemented `Reseeder` types

```@docs
NewLESReseeder
ReseedDisabled
```

