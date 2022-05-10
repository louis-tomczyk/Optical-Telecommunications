function out = AGD(in,learning_rate,memory_rate)

    assert(learning_rate>0 && learning_rate<=1,"\n\t learning rate must be in ]0,1]");
    assert(memory_rate>=0 && memory_rate<=1,"\n\t memory rate must be in [0,1]");
    
    alpha   = learning_rate;
    gamma   = memory_rate;

    if length(in) == 2
        vn = in(end)-in(end-1);
        gn = vn/alpha;
    else
        vn1     = in(end)-in(end-1);
        vn2     = in(end-1)-in(end-2);
        gn      = (vn1-gamma*vn2)/alpha;
    end

    out     = [in,in(end)+alpha*gn];
end