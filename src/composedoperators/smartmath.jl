"""
support that can be called on double or single composed integral operator. The function can make connections with the rest of BEAST, simplify the results etc
For now it just returns the original structure.
"""
smartmath(x) = x

smartmath(x::CompDoubleInt) = _CompDoubleInt(x.tfunc,t.pairing1,t.kernel,t.pairing2,t.bfunc)
smartmath(x::CompSingleInt) = _CompSingleInt(x.tfunc,t.pairing1,t.kernel,t.pairing2,t.bfunc)

#TODO no smartmath, remove this section we assume the user is "smart"