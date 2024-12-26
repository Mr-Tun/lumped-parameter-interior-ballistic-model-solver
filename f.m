function dydt = f(y)
global n rho_p omega f k alpha u1 e1 d chi lambda mu chi_s lambda_s m diameter V0 lg p0 phi theta delta S l0 rho ZK
    dydt = zeros(6, 1);
    dydt(2) = (u1/e1 * (y(5)^n)) * (y(2) < ZK);
    dydt(1) = (chi + 2*lambda*chi*y(2) + 3*mu*chi*y(2)^2) * dydt(2) * (y(2) < 1) + (chi_s/ZK * (1 + 2*lambda_s/ZK*y(2))) * dydt(2) * (y(2) >= 1 && y(2) < ZK);
    dydt(3) = y(4);
    dydt(4) = (S * y(5)) / (phi * m);
    dydt(6) = l0 * (-delta * (alpha - 1/rho_p) * dydt(1));
    dydt(5) = 1 / (S * (y(3) + y(6))) * (f * omega * dydt(1) - theta * phi * m * y(3) * dydt(4) - S * y(5) * (dydt(3) + dydt(6)));
end
