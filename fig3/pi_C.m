%% Purchase insurance income income
function E_C_two = pi_C(k,alpha,w,c,deta,UW,q,p,r,combMatrix,M,beta)
    if k>=M
        fuzhu1 = 0;
        fuzhu2 = 0;
        for s = 0:k
            for j = 0:k-s
                for m = 0:k-s-j
                    n = k - s - j - m;
                    if n==0
                        fuzhu1=fuzhu1+combMatrix(k+1,s+1)*combMatrix(k-s+1,j+1)*combMatrix(k-s-j+1,m+1)*((p-r)^s)*((1-q-r)^j)*((q+r-p)^m)*r^n*((s+j)*UW(w-beta*c-deta*w)+m*UW(w-beta*c+alpha*w-deta*w));
                        fuzhu2=fuzhu2+combMatrix(k+1,s+1)*combMatrix(k-s+1,j+1)*combMatrix(k-s-j+1,m+1)*((p-r)^s)*((1-q-r)^j)*((q+r-p)^m)*r^n*((s+j)*UW(w-beta*c)+m*UW(w-beta*c+alpha*w));
                    end
                end
            end
        end
        fuzhu = fuzhu2/k-fuzhu1/k;
        E_C_two=(1-q+p-2*r)*UW(w-beta*c-deta*w)+(q-p+r)*UW(w-beta*c+alpha*w-deta*w);
        for h=0:k-1
            E_C_two=E_C_two+r*nchoosek(k-1,h)*r^h*(1-r)^(k-1-h)*UW((1-alpha)*w-beta*c-deta*w+(k*w*deta)/(h+1));
        end
        E_C_two=E_C_two + fuzhu;
    else
        fuzhu1 = 0;
        fuzhu2=0;
        for s = 0:k
            for j = 0:k-s
                for m = 0:k-s-j
                    n = k - s - j - m;
                    if n==0
                        fuzhu1=fuzhu1+combMatrix(k+1,s+1)*combMatrix(k-s+1,j+1)*combMatrix(k-s-j+1,m+1)*((p-r)^s)*((1-q-r)^j)*((q+r-p)^m)*r^n*((s+j)*UW(w-c-deta*w)+m*UW(w-c+alpha*w-deta*w));
                        fuzhu2=fuzhu2+combMatrix(k+1,s+1)*combMatrix(k-s+1,j+1)*combMatrix(k-s-j+1,m+1)*((p-r)^s)*((1-q-r)^j)*((q+r-p)^m)*r^n*((s+j)*UW(w-c)+m*UW(w-c+alpha*w));
                    end
                end
            end
        end
        fuzhu = fuzhu2/k-fuzhu1/k;
        E_C_two=(1-q+p-2*r)*UW(w-c-deta*w)+(q-p+r)*UW(w-c+alpha*w-deta*w);
        for h=0:k-1
            E_C_two=E_C_two+r*nchoosek(k-1,h)*r^h*(1-r)^(k-1-h)*UW((1-alpha)*w-c-deta*w+(k*w*deta)/(h+1));
        end
        E_C_two=E_C_two + fuzhu;
    
    end

end
