module TrapRule
using LinearAlgebra

export traprule

function trap_step( f, x0, t0, t1; tol=1e-8, max_iter = 300)

    h = t1 - t0
    f0 = f(x0,t0)

    iter = 0
    flag = 0
    x1 = x0 + h*f0
    f1 = f(x1,t1)
    while flag == 0
        iter += 1

        x1 = x0 + h/2*( f1 + f0 )

        f1 = f(x1,t1)
        err = norm( - x1 + x0 + h/2*( f1 + f0 ) )

        if err ≤ tol
            flag = 1
        elseif iter ≥ max_iter
            flag = -1
        end
    end

    return x1, flag, iter
end

function traprule(f, t0, tf, h, x0)
    time = t0:h:tf
    nstep = length(time)
    nstates = length(x0)

    X = zeros(nstates, nstep)
    X[:,1] = x0
    for i = 2:nstep
        xi,flag,iter = trap_step(f,X[:,i-1],time[i-1],time[i])
        if flag ≠ 1
            println("Sub step failed to converge")
            break
        end
        X[:,i] = xi
    end

    return X,time
end


end
