function B0_unwrap = unwrap_B0_algorithms(TE,B0_wrapped,algoName)
switch algoName
    case 'Berglund'
        B0_unwrap=B0_unwrap_Fatty_Riot_GC(TE,B0_wrapped);
    case 'Snubben'
        B0_unwrap=B0_unwrap_Snubben(TE,B0_wrapped);
    case 'BerglundFWQPBO'
        B0_unwrap=B0_unwrap_BerglundFWQPB0(TE,B0_wrapped);
    otherwise
        warning('''algoName''Algorithm Name unknown')   
end
    
end


function B0_unwrapped=B0_unwrap_BerglundFWQPB0(TE,B0_wrapped)
        dt=diff(TE);
        dte=dt(1);
        B0_pi=pi*(B0_wrapped*2*dte-1);% [-pi, pi]
        B0_unwrapped=(unwrap3(B0_pi)/pi+1)/dte/2;
end


function B0_unwrapped=B0_unwrap_Snubben(TE,B0_wrapped)
        dt=diff(TE);
        dte=dt(1);
        B0_pi=pi*(B0_wrapped*2*dte-1);% [-pi, pi]
        B0_unwrapped=(unwrap3(B0_pi)/pi+1)/dte/2;
end


function B0_unwrapped=B0_unwrap_Fatty_Riot_GC(TE,B0_wrapped)
        dt=diff(TE);
        dte=dt(1);
        B0_pi=B0_wrapped*360*dte*(pi/180);% [-pi, pi]
        B0_unwrapped=unwrap3(B0_pi)*180/pi/dte/360;
end