function kappa = kappa4(vector)

y = 1-vector(1);
x = vector(2);
kappa =1;
m=1e4;

if abs(x-.5)<.4 && abs(y-.5)<.4
    if y > .1 && y < .2
        for i = 1:2:10
            if x > (i-1)*.1+.02 && x < (i-1)*.1+.07
                kappa = m*10;
            end
        end
    elseif y > .2 && y < .3
        for i = 1:2:10
            if x > (i-1)*.1+.04 && x < (i-1)*.1+.08
                kappa = m*10;
            end
        end
    elseif y > .4 && y < .5
        for i = 1:2:10
            if x > (i-1)*.1+.01 && x < (i-1)*.1+.05
                kappa = m*4;
            end
        end
    elseif y > .5 && y < .6
        for i = 1:2:10
            if x > (i-1)*.1+.02 && x < (i-1)*.1+.07
                kappa = m*4;
            end
        end
    elseif y > .6 && y < .7
        for i = 1:2:10
            if x > (i-1)*.1+.04 && x < (i-1)*.1+.08
                kappa = m*4;
            end
        end
    elseif y > .8 && y < .9
        for i = 1:2:10
            if x > (i-1)*.1+.01 && x < (i-1)*.1+.05
                kappa = m*6;
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end    
    if y > .2 && y < .3
        for i = 2:2:10
            if x > (i-1)*.1+.01 && x < (i-1)*.1+.05
                kappa = m*8;
            end
        end
    elseif y > .3 && y < .4
        for i = 2:2:10
            if x > (i-1)*.1+.02 && x < (i-1)*.1+.06
                kappa = m*8;
            end
        end
    elseif y > .4 && y < .5
        for i = 2:2:10
            if x > (i-1)*.1+.03 && x < (i-1)*.1+.07
                kappa = m*8;
            end
        end
    end
    for j = 1:10
        for i = 2:2:10
            if x > (i-1)*.1+.03 && x < (i-1)*.1+.06 && y > (j-1)*.1+.03 && y < (j-1)*.1+.06
                kappa = m*8;
            end
        end
    end
    if x > .4 && x < .6
        kappa = 1;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    f=inline('(0*sin(2*pi*x)+1.1*sqrt(x))*0.05+0.4');
    g=inline('(0*sin(2*pi*x)+1.1*sqrt(x))*0.05+0.42');
    
    if ceil(x*100)/100>=f(ceil(y*100)/100) && ceil(x*100)/100<=g(ceil(y*100)/100)
        kappa = m*10;
    end
end

