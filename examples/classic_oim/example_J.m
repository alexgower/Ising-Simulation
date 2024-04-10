function J = example_J()

    J = zeros(8,8);

    % MAXCUT on 8 Nodes J_ij (Listing 1.2) ---
    J(1, 2) = -1;
    J(2, 3) = -1;
    J(3, 4) = -1;
    J(4, 5) = -1;
    J(5, 6) = -1;
    J(6, 7) = -1;
    J(7, 8) = -1;
    J(1, 8) = -1;
    J(1, 5) = -1;
    J(2, 6) = -1;
    J(3, 7) = -1;
    J(4, 8) = -1;

    J = J + transpose(J);

end