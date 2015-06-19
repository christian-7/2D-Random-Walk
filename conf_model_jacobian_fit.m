function  [fitParam,resnorm,residual,exitflag]=conf_model_jacobian_fit(par0,lowerBound, upperBound,x,y)

    options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on','Display','off');
  
    [fitParam,resnorm,residual,exitflag] = lsqnonlin(@(p) ...
            f_confined_model(p, x,y),...
            par0, lowerBound, upperBound, options);
end

 function [err, jacobian] = f_confined_model(p,x,y)

        fun = p(1).^2*(1-p(2)*exp((-4.*p(3).*p(4).*x)/p(1).^2));
      
        err = reshape(fun-y,1,[]);
  
   if nargout > 1
       jacobian = f_confined_model_Jacobian(p, x);
   end
   
 end
 
 function jacobian = f_confined_model_Jacobian(p, x)
% the Jacobian J is an m-by-n matrix where J(i,j) is the partial derivative 
% of F(i) with respect to par(j). (The Jacobian J is the transpose of the gradient of F.)
% F is the ouput of the function Zmodel that accepts a vector of parameters par. 

jacobian1 = 2.*p(1) ...
            -2.*p(1).*p(2).*exp((-4.*p(3).*p(4).*x)/(p(1).^2)) ...
            -p(1).^2*p(2).*exp((-4.*p(3).*p(4).*x)/(p(1).^2)).*((8.*p(3).*p(4).*x)/(p(1).^3));
       
jacobian2 = -p(1).^2*exp((-4.*p(3).*p(4).*x)/(p(1).^2));

jacobian3 = -p(1).^2*p(2).*exp((-4.*p(3).*p(4).*x)/(p(1).^2)).*((-4.*p(4).*x)/(p(1).^2));
jacobian4 = -p(1).^2*p(2).*exp((-4.*p(3).*p(4).*x)/(p(1).^2)).*((-4.*p(3).*x)/(p(1).^2));

jacobian = [reshape(jacobian1,[],1) reshape(jacobian2,[],1) ...
          reshape(jacobian3,[],1) reshape(jacobian4,[],1)];
          
 end