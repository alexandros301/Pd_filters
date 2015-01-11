filters is a collection of Pure Data external filter objects based on Mike Moser-Booth's [filtercoeff.mmb~] and [biquad.mmb~] abstractions, translated in C.
Instead of the [biquad.mmb~] abstraction, the proposed use of [fexpr~] is used in the code.

The library includes:
[lowPass~]
[highPass~]
[bandPass~]
[allPass~]
[resonant~]
[lowShelf~]
[highShelf~]
[peakNotch~]

Help patches are included. The help patches use the [spectrum.mmb~] abstraction for spectrum projection, as this is used in [filtercoeff.mmb~]'s 
abstraction help patch.

Written by Alexandros Drymonits
alexandros[at]drymonitis.me

July 2014
