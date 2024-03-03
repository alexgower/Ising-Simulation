function [Ebest,flag] = breakout(E,Ebest,Esol,record)

flag = 0;


if E > Ebest
    % Update Ebest if E is better (remember Pei uses higher energy better convention)
    Ebest = E;

    % Break out of the loop if there exists a solution energy (Esol !=0), and we're not in record mode
    % and the difference between the best energy and the solution energy is less than 0.01
    if abs(Esol-Ebest)<0.01 && ~record && Esol
        flag = 1;
    end
end

end