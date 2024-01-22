# Notes on how to make the case ready for simulation

**After initialization the case using `prepareBeams` (it will initialize the beam posture and copy the data generated to CFD case), make sure to apply the followings manually:**

### 1- changing W boundary condition

For beamone, the `W` (e.g. for **beamone**, it is located in: *`mooredFloatingObject/0.orig/beamone/W`* ) looks like this:

```cpp
boundaryField
{
    left
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
    right
    {
        displacementSeries
        {
            file            "$FOAM_CASE/constant/timeVsDisplacement";
            outOfBounds     clamp;
        }
        type            fixedDisplacement;
        value           uniform (-0.16093 -0.120518 -0.100692);
    }
    beam
    {
        type            empty;
    }
}
```
The `right` BC for **beamone** and **beamthree** needs to be changed to `fixedValue` like below:
```cpp
boundaryField
{
    left
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
    right
    {
        type            fixedValue;
        value           uniform (-0.16093 -0.120518 -0.100692);
    }
    beam
    {
        type            empty;
    }
}
```

Since in **beamtwo** and **beamfour**  both `left` and `right` patches are using `timeVsDisplacement`, make sure to apply similar modification to both `left` and `right` patches.

___

### 2- changing Theta boundary condition
The path to `timeVsMoment` needs to be updated, the `Theta` (for **beamone**, it is located in: located in *`mooredFloatingObject/0.orig/beamone/Theta`* ) looks like this:

```cpp
boundaryField
{
    left
    {
        type            momentBeamRotationNR;
        value           uniform (-0.18179 0.430172 -0.0480808);
        momentSeries
{
        file            "$FOAM_CASE/constant/timeVsMoment";
        outOfBounds     clamp;
}
        moment          uniform (0 0 0);
    }
    right
    {
        type            momentBeamRotationNR;
        value           uniform (-0.254974 -0.885381 -0.170938);
        momentSeries
{
        file            "$FOAM_CASE/constant/timeVsMoment";
        outOfBounds     clamp;
}
        moment          uniform (0 0 0);
    }
    beam
    {
        type            empty;
    }
}
```

which the path to *`timeVsMoment`* should be updated like (adding regions middle directory):

```cpp
boundaryField
{
    left
    {
        type            momentBeamRotationNR;
        value           uniform (-0.18179 0.430172 -0.0480808);
        momentSeries
{
        file            "$FOAM_CASE/constant/beamone/timeVsMoment";
        outOfBounds     clamp;
}
        moment          uniform (0 0 0);
    }
    right
    {
        type            momentBeamRotationNR;
        value           uniform (-0.254974 -0.885381 -0.170938);
        momentSeries
{
        file            "$FOAM_CASE/constant/beamone/timeVsMoment";
        outOfBounds     clamp;
}
        moment          uniform (0 0 0);
    }
    beam
    {
        type            empty;
    }
}
```

Make sure to update this path for each beam in the solution. (e.g. *`"$FOAM_CASE/constant/beamtwo/timeVsMoment";`* for **beamtwo**)


Next step is to execute the `Allrun` script located in *`mooredFloatingObject/`* and wait for the results.

____

**Note**: The `cleanBeams` script will remove beam meshes in *`cable_initialization`* and beam regions in *`mooredFloatingObject`*.