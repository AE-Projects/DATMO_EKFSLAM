function plot_SLAM_DATMO(map, true_pose,z_t_s, z_t_m, true_mov, cov_matrix_all, mu_all, cov_mov, mu_mov, color)
[two, I] = size(z_t_s);
cla;
scatter(map(1,:), map(2,:))
hold on
grid on
mu_j_t = zeros(2,I);
[two, mov_sz] = size(z_t_m);
if mov_sz > 0
    mu_mov2 = measurement_model_inverse(z_t_m, true_pose); 
    plot([true_pose(1,1), mu_mov2(1,1)],[true_pose(2,1), mu_mov2(2,1)], 'r-','LineWidth',0.3)
end

for i = 1:I
    mu_j_t(:,i) = measurement_model_inverse(z_t_s(:,i), true_pose);    
end
plot(true_pose(1,1), true_pose(2,1), 'k+', 'MarkerSize', 5, 'LineWidth', 1.5);
plot(true_mov(1,:), true_mov(2,:), 'ro', 'MarkerSize', 5, 'LineWidth', 1.5);
plot([true_pose(1,1)*ones(1,I);[mu_j_t(1,:)]],[true_pose(2,1)*ones(1,I);[mu_j_t(2,:)]],'g-','LineWidth',0.3);
%pause(0.5)

[n_ellipse, two] = size(mu_all);

for n = 1:n_ellipse+1 
    
    if n==n_ellipse+1
        if mu_mov==0
            break
        end
        cov_matrix = cov_mov(1:2, 1:2);
        origin = mu_mov(1:2,1)';
    else
        cov_matrix = cov_matrix_all(2*n-1:2*n, 2*n-1:2*n);
        origin = mu_all(n,1:2);
        cov_matrix = 100*cov_matrix;
    end
    
    [eigT, eigV] = eig(cov_matrix*cov_matrix'); 
    eigenvalues = diag(eigV);
    if any(eigenvalues <= 0)
        error('Covariance matrix must be positive definite.');
    end
    a = sqrt(eigenvalues(1)); % Semi-major axis
    b = sqrt(eigenvalues(2)); % Semi-minor axis
    theta = linspace(0, 2*pi, 100);
    ellipse_points = eigT * [a 0; 0 b] * [cos(theta); sin(theta)];
    ellipse_points(1, :) = ellipse_points(1, :) + origin(1);
    ellipse_points(2, :) = ellipse_points(2, :) + origin(2);
    if n==1
        color_ = 'g';
    elseif n==n_ellipse+1
        color_ = 'cyan';
    else
        color_ = color;
    end
    plot(ellipse_points(1, :), ellipse_points(2, :), 'Color', color_, 'LineWidth', 1.5);
    %hold on;
    plot(origin(1), origin(2), 'k+', 'MarkerSize', 5, 'LineWidth', 1.5); 
    axis equal; 
    %axis([-1 12 -1 12])
    grid on;
    
end


end
