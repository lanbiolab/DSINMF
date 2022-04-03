# CaFew_Project
cluster-aware feature weighting
--------------------------------
This project includes 1 main function and 7 subfunctions
1. InitCenter.r : The function choose the initial cluster center for data;
2. init.r: The function returns the initialization value of parameter delte and weight matrix initially;
3. UpdateWeight.r: The function returns the value of updated weight matrix in certain iteration;
4. Patition.r: The function return the patition results of samples in certain iteration;
5. UpdateCenter.r: The function returens the updated cluster centers in certain iteration;
6. UpdateDelta.r: The function returens the updated parameter delta in certain iteration;
7. objFun.r: The function returns the value of objective function in certain iteration;
8. Main.r: The mian function involk above 7 subfunctions and selected cluster-aware features for data.
---------------------------------
To use this method, please read your data format as cells in rows and genes in colums and then run the file Main.r
