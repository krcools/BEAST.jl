struct HH3DFarField{K}
    gamma::K
end
  
function HH3DFarField(;wavenumber=error("wavenumber is a required argument"))
    if iszero(real(wavenumber))
        HH3DFarField(-imag(wavenumber))
    else
        HH3DFarField(wavenumber*im)
    end
end
  
quaddata(op::HH3DFarField,rs,els) = quadpoints(rs,els,(3,))
quadrule(op::HH3DFarField,refspace,p,y,q,el,qdata) = qdata[1,q]

kernelvals(op::HH3DFarField,y,p) = exp(op.gamma*dot(y,cartesian(p)))
function integrand(op::HH3DFarField,krn,y,f,p)
    krn * f.value
end