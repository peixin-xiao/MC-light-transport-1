this folder fix some problems in previous code : 
in previous code, every time after scatter of refraction, the coordinate x,y,z will :
  x = x + length_voxel * n
  this is not precise, instead of it, it should be : 
  x = x + step * 1
  (step = time_delta * LIGHT_SPEED
when the scattering coefficient is very high, the result seems unresonable and time consuming
