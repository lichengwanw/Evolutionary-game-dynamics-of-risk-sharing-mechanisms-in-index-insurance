%% Purchase insurance fitness calculation
%%The results are consistent with those calculated using p(u,v,s,t) and q(u,v,s,t).
function F_C_CII_last = F_C_CII(i,Z,N,alpha,w,c,deta,UW,q,p,r,combMatrix,M,beta)
   
    function E_C_two = E_C_CII(k) 
        if k>=M
            E_C_two=(1-q+p-2*r)*UW(w-beta*c-deta*w)+(q-p+r)*UW(w-beta*c+alpha*w-deta*w);
            for h=0:k-1
                E_C_two=E_C_two+r*combMatrix(k,h+1)*r^h*(1-r)^(k-1-h)*UW((1-alpha)*w-beta*c-deta*w+(k*w*deta)/(h+1));
            end
        else
            E_C_two=(1-q+p-2*r)*UW(w-c-deta*w)+(q-p+r)*UW(w-c+alpha*w-deta*w);

            for h=0:k-1
                E_C_two=E_C_two+r*combMatrix(k,h+1)*r^h*(1-r)^(k-1-h)*UW((1-alpha)*w-c-deta*w+(k*w*deta)/(h+1));
            end

        end

    end
    
    F_C_CII_last=0;
    for k=0:N-1
            if i-1<k || (Z-i)<(N-1-k)
                F_C_CII_last=F_C_CII_last+0;
            else

                F_C_CII_last=F_C_CII_last+E_C_CII(k+1)*combMatrix(i,k+1)*combMatrix(Z-i+1,N-k)/combMatrix(Z,N);
            end
            
        
    end
end