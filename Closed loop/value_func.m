%% Value function
function V = value_func(beta_plus,beta_minus,lambda,u,R)
    
    gt = u >= R;
    V = zeros(size(u));
    V(gt) = (u(gt)-R(gt)).^beta_plus(gt);
    V(~gt) = -lambda(~gt).*((R(~gt)-u(~gt)).^beta_minus(~gt));
    
%     if (u >= R)
%         V = (u-R).^beta_plus;
%     else
%         V = -lambda.*((R-u).^beta_minus);
%     end

end