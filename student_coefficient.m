function t = student_coefficient(NPoints,ConfidenceRate,Coeffs)

    assert(NPoints>=3 && NPoints<=120,"\n\t Please take at least 3 measures and maxmimum")
    assert(ConfidenceRate==95||ConfidenceRate==99,"\n\t Available confidence rates are 95% or 99%")

    switch ConfidenceRate
        case 95
            if NPoints < 120
                t = Coeffs(NPoints-2,2);
            else
                t = Coeffs(end,2);
            endif
        case 99
            if NPoints < 120
                t = Coeffs(NPoints-2,3);
            else
                t = Coeffs(end,3);
            endif
    endswitch
endfunction