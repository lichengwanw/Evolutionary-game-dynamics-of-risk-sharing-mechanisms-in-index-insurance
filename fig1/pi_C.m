%% Purchase insurance income income
function E_C_two = pi_C(k,alpha,w,c,deta,UW,q,p,r,combMatrix,M,beta)
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
