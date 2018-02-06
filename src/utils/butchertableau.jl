using StaticArrays



"""
    butcher_tableau_radau_2stages()

Returns (A,b,c) corresponding to the Butcher tableau for the 2 stage Radau IIA scheme.
"""
function butcher_tableau_radau_2stages()
	A = @SMatrix [5/12   -1/12
	              3/4    1/4];
	b = @SVector [3/4, 1/4];
	c = @SVector [1/3, 1.0];
	return (A, b, c);
end

"""
    butcher_tableau_radau_3stages()

Returns (A,b,c) corresponding to the Butcher tableau for the 3 stage Radau IIA scheme.
"""
function butcher_tableau_radau_3stages()
	A = @SMatrix [(88.0-7.0*sqrt(6.0))/360.0       (296.0-169.0*sqrt(6.0))/1800.0   (-2.0+3.0*sqrt(6.0))/225.0
	              (296.0+169.0*sqrt(6.0))/1800.0   (88.0+7.0*sqrt(6.0))/360.0       (-2.0-3.0*sqrt(6.0))/225.0
	              (16.0-sqrt(6.0))/36.0            (16.0+sqrt(6.0))/36.0             1/9];
	b = @SVector [(16.0-sqrt(6.0))/36.0,    (16.0+sqrt(6.0))/36.0,    1/9];
	c = @SVector [(4.0-sqrt(6.0))/10.0,     (4.0+sqrt(6.0))/10.0,     1.0];
	return (A, b, c);
end
