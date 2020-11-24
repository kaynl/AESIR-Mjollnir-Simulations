function [output] = Wall_ext_temp_finder(T_ext_wall,T_boundary,h,opts)
%EXTERNAL WALL TEMPERATURE FINDER by getting the zero of this function
%http://dark.dk/documents/technical_notes/simplified%20aerodynamic%20heating%20of%20rockets.pdf

sigma = opts.stephan_cst;
epsilon = opts.aluminium_emissivity;


output = h - (sigma*epsilon*T_ext_wall^4)/(T_boundary-T_ext_wall);

end

