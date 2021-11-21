function [obj,x,status,loop] = OneNormLP3037534676(A,b)
%This function solves a one_norm regression problem
%obj denotes the objective value 
%x denotes the x s.t |Ax-b| achieve the minimum value
%status contains several messages:1)degenerate
%2)optimal reached --> obj converges
%There will be no other conditions for the status since we only have 4
%situations from duality theorem and the situation of the primal is
%feasible but without a lower bound is removed because of one_norm. It's
%also feasible because we don't impose any restriction on x. Thus, it must
%be feasible. The only 2 conditions are whether it's degenerate. 

 %setup initial value
    [m,n] = size(A);
    obj = -1;
    x = NaN*ones(n,1);
    status = 'initial state';
    loop = 0;
    opt_flg = 0;
    %we start with a random basis
    %we add phase I to deduce the possibility of failing 
    %B is a row vector
    
    %remove this for large case. 
    %{
    all_comb = nchoosek(1:1:m,n);
    for i=1:nchoosek(m,n)
        B = all_comb(i,:);
        if det(A(B,:))==0
            status = 'degenerate';
            return;
        else
            x_buffer = inv(A(B,:))*b(B,:);
            check = A*x_buffer-b;
            check_row = find(check<1e-10);
            if size(find(check(check_row,:)>-1e-10),1)>n
                status = 'degenerate';
                return;
            end
        end
    end
    %}
    %    
    k = randperm(m);
    B = k(1:n);
    M = inv(A(B,:));
    h = zeros(m,1);
    y = zeros(m,1);
    t = zeros(m,1);
    obj_buffer = 0;
    if det(A(B,:))==0
        status = 'degenerate';
        return;
    end
    while opt_flg == 0
        B_bar = setdiff(k,B);
        % This determine whether the A(B,:) is singular
        if det(A(B,:)) == 0 || cond(A(B,:))>1e3
            status = 'degenerate';
            return
        end
        x = M*b(B,:);
        h=A*x-b;
        % the obj is supposed to be less than the previous one based on
        % non-degeneracy assumption. If it voilates the condition, it's
        % degenerate
        if loop == 0
            obj_buffer = sum(abs(h));
        else
            if obj_buffer < sum(abs(h))
                status = 'degenerate';
                x = NaN;
                return
            end
        end
        loop = loop + 1;
        disp("The new basis is")
        disp(B);
        disp("The current objective value is " + sum(abs(A*x-b)));
        if size(find(abs(h(B_bar,:)) < 1e-10),1) > 0 
            status = 'degenerate';
            return
        end
        y(B_bar,:)=sign(h(B_bar,:));
        y(B,:) = -transpose(M)*(transpose(A(B_bar,:))*y(B_bar,:));
        y_B = y(B,:);
        y_abs = abs(y(B,:));
        if size(y_B(y_abs>1),1) == 0
            opt_flg = 1;
            status = 'done!';
            obj = sum(abs(A*x-b));
            return;
        end
        %choose s
        %y_s_choice = y_B(y_abs>1);
        row_choice = find(y_abs>1);
        s = row_choice(1);
        y_s = y_B(s);
        %choose  r
        t(B_bar,:) = -(sign(y_s))*(y(B_bar,:)).*(A(B_bar,:)*M(:,s));
        t_bar = t(B_bar,:);
        h_bar = h(B_bar,:);
        h_j = h_bar(t_bar>0);
        t_j = t_bar(t_bar>0);
        h_divide_t = abs(h_j)./t_j;
        [~,row] = min(h_divide_t);
        real_tj = t_j(row);
        r = B_bar(t_bar==real_tj);
        B(B==B(s)) = r;
        % Sherman-Morrison
        theta = transpose(A(r,:)*M);
        M(:,s)=(1/(theta(s)))*M(:,s);
        for i=1:n
            if i~= s
                M(:,i)=M(:,i)-theta(i)*M(:,s);
            end
        end
        
        
    end