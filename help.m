for iteration = 1:numberofit
    clear densfieldnew
    densfieldnew = zeros(1,cellsperdim);
    
    for rad = 1:cellsperdim
        clear denshere densout
        denshere = densfield(rad);
        if particlefield(rad)
            densout = denshere - particleparts;
            densfieldnew(rad) = densfieldnew(rad) + particleparts; % for both (1) and (2)!
        else, densout = denshere;
        end
        if densout <= 0, continue; end
        
        for constituent = 1:densout
            clear velo dist_rad dist time
            velo = initvelo*cellsperdim*(1-min(denshere/densmax,1));
            dist_rad = cos(rand*pi);
            dist = 0; % dist = dist_rad;
            time = 0;
            
            while time < 1 && velo > 0 % if velo <= 0, continue, so don't enter BH
                dist = dist + dist_rad;
                time = time + (1/velo); % t' = t + dt = t + dx/v
                clear new_rad velo
                new_rad = round(mod(rad + dist,cellsperdim)); if new_rad == 0, new_rad = cellsperdim; end
                velo = initvelo*cellsperdim*(1-min(densfield(new_rad)/densmax,1));
            end
            
            % add radial acceleration to particle location
            new_rad = round(mod(rad + dist + accelerate,cellsperdim)); if new_rad == 0, new_rad = cellsperdim; end
            
            if densfieldnew(new_rad) < densmax
                densfieldnew(new_rad) = densfieldnew(new_rad)+1;
            elseif densfieldnew(rad) < densmax % stay put if you intend to add to new BH area
                densfieldnew(rad) = densfieldnew(rad)+1;
            end
        end
    end
    % if any(densfieldnew > densmax), return; end
    
    clear densfield
    if gradmin < gradmax
        densfield = densgrad;
    else
        densfield = densfieldnew;
    end
    % linear rescale
    particle = densfieldnew(particlefield);
    particlegrad = particle(end)-particle(1);
    particle_new = particle(1) + (0:(length(particle)-1))*particlegrad/(length(particle)-1);
    densfield(particlefield) = round(particle_new/sum(particle_new)*(n_tot - sum(densfield(~particlefield))));
    
    % averaging
    densfield_tot = densfield_tot + densfield;
end