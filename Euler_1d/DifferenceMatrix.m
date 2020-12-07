function DM = DifferenceMatrix(datacoord, centercoord)

[dr, cc] = ndgrid(datacoord(:), centercoord(:));

DM = dr - cc;