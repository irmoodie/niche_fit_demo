data{
    int N;
    vector[10] r0_sd;
    vector[10] r0;
    int temperature[10];
    vector[10] tref;
    vector[10] k;
}
parameters{
    vector[N] r0_true;
    real rtref;
    real th;
    real e;
    real eh;
    real<lower=0> sigma;
}
model{
    vector[10] est;
    sigma ~ exponential( 1 );
    eh ~ normal( log(3) , 0.4 );
    e ~ normal( log(0.5) , 0.3 );
    th ~ normal( 35 , 3 );
    rtref ~ normal( log(0.2) , 0.3 );
    for ( i in 1:10 ) {
        est[i] = (exp(rtref) * exp(exp(e)/k[i] * (1/tref[i] - 1/(temperature[i] + 273.15)))) * (1/(1 + exp(exp(eh)/k[i] * (1/(th + 273.15) - 1/(temperature[i] + 273.15)))));
    }
    r0_true ~ normal( est , sigma );
    r0 ~ normal( r0_true , r0_sd );
}