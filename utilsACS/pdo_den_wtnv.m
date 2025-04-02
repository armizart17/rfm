function [x, cost, error, fid, reg] = pdo_den_wtnv(y, lambda, tau, maxIter, tol, stableIter, weightEstimators)
% function [x, cost, error, fid, reg] = pdo_den_wtnv(y, lambda, tau, maxIter, tol, stableIter, weightEstimators)
	
	rho = 1.99;		% relaxation parameter, in [1,2]
	sigma = 1/tau/8; % proximal parameter
	[H,W,C]=size(y);
    
    global weights
    weights = weightEstimators;
 
 	x2 = y; 		% Initialization of the solution
       
	u2 = zeros([size(y) 2]); % Initialization of the dual solution
	cy = sum(sum(sum(y.^2)))/2;
	primalcostlowerbound = 0;
		
	ee = 1;	
    error(1) = 1;
    iter = 1;
    cost(1) = 100;
    fid(1) = 1;
    reg(1) = 1;
    
    disp('Itn   |    Cost   |   Delta(Cost)');
    disp('----------------------------------');
    
    while (iter < maxIter) && (ee > tol)  
		x = prox_tau_f(x2-tau*opDadj(u2),y,tau);
		u = prox_sigma_g_conj(u2+sigma*opD(2*x-x2), lambda);
        
        
		x2 = x2+rho*(x-x2);
		u2 = u2+rho*(u-u2);
        
        fid(iter+1) = 0.5*sum(sum(sum((x-y).^2)));
        reg(iter+1) = lambda*TVnuclear_EMZ(x);

%         cost(iter+1) = 0.5*sum(sum(sum((x-y).^2)))+lambda*TVnuclear_EMZ(x);
        cost(iter+1) = fid(iter+1) + reg(iter+1);

%         cost(iter+1) =  abs ( cy-sum(sum(sum((y-opDadj(u)).^2)))/2 ); % 2nd order 

        ee = abs(cost(iter+1) - cost(iter));     
        % NO NORMALIZADO
%         error(iter+1) = ee / cost(iter+1);
        
        % NORMALIZADO
        ee = ee / cost(iter+1);
        error(iter+1) = ee;
        
        if (iter < stableIter) % en las primeras 100 iteraciones inestable 
            ee = 1; % poner error grande (cualquier valor) para que no pare iteracion            
        end
        
		if mod(iter,5)==0
			primalcost = 0.5*sum(sum(sum((x-y).^2)))+lambda*TVnuclear_EMZ(x);
            
			dualcost = cy-sum(sum(sum((y-opDadj(u)).^2)))/2;

			primalcostlowerbound = max(primalcostlowerbound,dualcost);

% 			fprintf('%4d  | %f | %f | %e\n',iter,primalcost,...
% 				primalcostlowerbound,primalcost-primalcostlowerbound);
%             
            % fprintf('%4d  | %f | %e\n',iter, cost(iter+1), error(iter+1));
            fprintf('%4d  | %.3f | %.3fe-3\n', iter, cost(iter+1), error(iter+1)/1e-3);

%%%% For time execution

%             reshape_x = reshape(x, [], C);
%             snr_ratios = abs(mean(reshape_x, 1, 'omitnan')) ./ std(reshape_x, 0, 1, 'omitnan');
            %%%% For time execution
% 
%             figure(112), 
%             plot(snr_ratios), title('SpRatios SNR'), grid;
%             ylim([0 4]), yticks(0:1:4)

        end
        iter = iter+1;
	end
return

function [u] = prox_tau_f(v, B, tau)
    
    u =  (v+tau*B)/(1+tau);

return

function [u] = prox_sigma_g_conj(v, lambda)

    u = v - prox_nuc_weight(v,lambda);
return

function [u] = opD(v)
    % Description: Dy(:,:,1) Dx (:,:,:,2), 3rd dimension is channel
    [H,W,C]=size(v);
    % usually u = [Dx(v), Dy(v)]; for each channel use for
    
%     jacobian_vImage = DxDy(vectorial_image(:,:,1));
% 
%     for i = 2:size(vectorial_image,3)
%         v = vectorial_image(:,:,i);
%         jacobian_vImage = [jacobian_vImage; DxDy(v)];
%     end
    
    % one line
%    u = cat(4,[diff(v,1,1);zeros(1,W,3)],[diff(v,1,2) zeros(H,1,3)]);
    u = cat(4,[diff(v,1,1);zeros(1,W,C)],[diff(v,1,2) zeros(H,1,C)]);
return

function [u] = opD_adap(v, factor)
    % Description: Dy(:,:,1) Dx (:,:,:,2), 3rd dimension is channel
    [H,W,C]=size(v);
    % usually u = [Dx(v), Dy(v)]; for each channel use for
    
%     jacobian_vImage = DxDy(vectorial_image(:,:,1));
% 
%     for i = 2:size(vectorial_image,3)
%         v = vectorial_image(:,:,i);
%         jacobian_vImage = [jacobian_vImage; DxDy(v)];
%     end
    
    % one line
%    u = cat(4,[diff(v,1,1);zeros(1,W,3)],[diff(v,1,2) zeros(H,1,3)]);
%     u = cat(4,[diff(v,1,1);zeros(1,W,C)],[diff(v,1,2) zeros(H,1,C)]);
    
    izq = [diff(v,1,1);zeros(1,W,C)] .* factor;
    der = [diff(v,1,2) zeros(H,1,C)] .* factor;
    u = cat(4, izq, der);
    
return



function [u] = opDadj(v)
    % u = Dy Dx
    u = -[v(1,:,:,1);diff(v(:,:,:,1),1,1)]-[v(:,1,:,2) diff(v(:,:,:,2),1,2)];	
return
    
function [u] = opDadj_adap(v, factor) %% ADAPTATIVO
    % u = Dy Dx
%     u = -[v(1,:,:,1);diff(v(:,:,:,1),1,1)]-[v(:,1,:,2) diff(v(:,:,:,2),1,2)];
    
    
    u = factor .* (-[v(1,:,:,1);diff(v(:,:,:,1),1,1)]- ...
        [v(:,1,:,2) diff(v(:,:,:,2),1,2)]);
return


function [u] = TVnuclear_EMZ(v)

    % u = norm(eig(opD(v)),1); % does not work 
    Jacobian = opD(v);
    u = NuclearNorm(Jacobian); %||J(u)||_nuclear
    
return


%%%% internet function %%%%%%%% nuclear matrix of matrix N dim
function val = NuclearNorm(y)

	s = diff(sum(y.^2,3),1,4);
	theta = atan2( 2*dot(y(:,:,:,1), y(:,:,:,2),3), -s )/2;
	c = cos(theta);
	s = sin(theta);
    
	val = sum(sum(sqrt(sum((bsxfun(@times,y(:,:,:,1),c)+... 
		bsxfun(@times,y(:,:,:,2),s)).^2,3)),2),1)+...
		sum(sum(sqrt(sum((bsxfun(@times,y(:,:,:,2),c)-...
		bsxfun(@times,y(:,:,:,1),s)).^2,3)),2),1);
    
    % EQUIVALENTE
%     val2 = sum( sum(sqrt( sum( ( y(:,:,:,1).*c + y(:,:,:,2).*s ).^2, 3) ), 2), 1 )+...
% 		   sum( sum(sqrt( sum( ( y(:,:,:,2).*c - y(:,:,:,1).*s ).^2, 3) ), 2), 1 );
   
return

function xx = prox_nuc_weight(y, lambda)
    eee = 0.00001;
    global weights;
    % GIVENS ROTATION (instead of SVD decomposition)
	s = diff( sum(y.^2,3), 1, 4 );
    
	theta = atan2( 2*dot(y(:,:,:,1), y(:,:,:,2),3), -s )/2;
    
	c = cos(theta);
    s = sin(theta);
    
    % S_matrix = [c -s; 
    %             s  c];
    
    x = cat( 4, y(:,:,:,1).*c + y(:,:,:,2).*s, ...
		       -y(:,:,:,1).*s + y(:,:,:,2).*c ); % x diagonalization 
    
    for ii = 1 : size(y,3)         
%         x(:,:,ii,:) = x(:,:,ii,:)* ((weights(ii))); 
       x(:,:,ii,:) = x(:,:,ii,:)* (1./ (weights(ii)+eee));
    end
           
    % PROXIMAL OPERATOR l2 ball
	tmp = max( sqrt(sum(x.^2,3)), lambda );
    
    tmp =  x .* ( (tmp-lambda)./tmp );
    
    % S_matrix_T = [c  s;
    %               -s c]
    %     xx = tmp*S_matrix_T 
    
    % RETURN FINAL RESULT dot product by vectors dim1 
	xx = cat(4, tmp(:,:,:,1).*c - tmp(:,:,:,2).*s, ...
		        tmp(:,:,:,1).*s + tmp(:,:,:,2).*c );
return

function x = prox_g(y, lambda)

	    % GIVENS ROTATION (instead of SVD decomposition)
	s = diff( sum(y.^2,3), 1, 4 );
    
	theta = atan2( 2*dot(y(:,:,:,1), y(:,:,:,2),3), -s )/2;

	c = cos(theta);
	s = sin(theta);

    % S_matrix = [c -s; 
    %             s  c];

	x = cat( 4, y(:,:,:,1).*c + y(:,:,:,2).*s, ...
		       -y(:,:,:,1).*s + y(:,:,:,2).*c ); % x diagonalization 
    
	% PROXIMAL OPERATOR l2 ball
	tmp = max( sqrt(sum(x.^2,3)), lambda );

	tmp =  x .* ( (tmp-lambda)./tmp );

    % S_matrix_T = [c  s;
    %               -s c]
    %     xx = tmp*S_matrix_T 

	x = cat(4, tmp(:,:,:,1).*c - tmp(:,:,:,2).*s, ...
		       tmp(:,:,:,1).*s + tmp(:,:,:,2).*c );
return