Note that angels shoud be like this (otherwise there is no garantee the code works, actually it does not work):

-90  < DEC < 90
-180 <  RA < 180

-> 

phi   = RA
theta = 90 - DEC

->

long = RA
lati = DEC

this structure is CONSITENT with HEALPY and it works verywell.
