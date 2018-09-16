% ransac show figure
n = length(A); 
A2 = (R * A') + repmat(t, 1, n); 
A2 = A2'; 
err = A2 - B;
err = err .* err;
err = sum(err(:));
rmse = sqrt(err/n)
% disp(sprintf("RMSE: %f", rmse));
% disp("If RMSE is near zero, the function is correct!");
% figure(1);clf 
% pcshow(A2, 'r'); hold on; pcshow(B,'b');