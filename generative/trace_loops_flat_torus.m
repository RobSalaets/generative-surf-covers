function [loop, loopps] = trace_loops_flat_torus(fl, gl, lm_idx, dir, prev_loop)
    %Find loops in plane space, use connectivity from uncut torus, lm_idx
    %in (cut) plane space, prev_loop in uncut torus space
    A = adjacency_matrix(gl.T_torus);
    if ~isempty(prev_loop)
        prev_loop = setdiff(prev_loop, fl.I_cut_to_uncut(lm_idx));
    end
    V = fl.V_plane;
    loop(1) = lm_idx;
    i = 1;
    curr = lm_idx;
    while curr ~= lm_idx || length(loop)  == 1
        Acurr = fl.I_cut_to_uncut(curr);
        [~, Aring] = find(A(Acurr,:));
        if ~isempty(prev_loop)
            Aring = setdiff(Aring, prev_loop);
        end
        if length(loop) > 1
            Aprev = fl.I_cut_to_uncut(loop(end-1));
            Aring = setdiff(Aring, Aprev);
        end
        if length(loop) == 200
            disp prijs
        end
        [ring, aj] = find(fl.I_cut_to_uncut == Aring);
        %Deal with boundary
        ringorig = ring;
        if length(ring) > length(Aring)
            counts = sum(aj==aj');
            for ii=1:length(counts)
                if counts(ii) > 1
                    all_ids = find(aj==aj(ii));
                    [~, mloc] = min(V(ring(all_ids),dir));
                    ring(all_ids) = ring(all_ids(mloc)); %edit ring to contain min along dir
                end
            end
        end
        if dir==1
            vs = mod(V(ringorig,:)-repmat([V(curr,1)-0.5 V(lm_idx,2)-0.5], [length(ringorig) 1]),1)...
                 + repmat([-0.5 -0.5], [length(ringorig) 1]);
        else
            vs = mod(V(ringorig,:)-repmat([V(lm_idx,1)-0.5 V(curr,2)-0.5], [length(ringorig) 1]),1)...
             + repmat([-0.5 -0.5], [length(ringorig) 1]);
        end
        vn = vecnorm(vs')';
        vs = vs./ vn;
        cost = (vs(:,dir) - 1).^2 + 4*vs(:,3-dir).^2;
        [~,iloc] = min(cost);
        i = i+1;
        loop(i) = ring(iloc);
        curr = loop(i);
        if length(loop) > size(gl.V_torus, 1)
            assert(0, "Loop cutting failed")
        end
    end
    loopps = loop;
    loop = fl.I_cut_to_uncut(loop);
end
