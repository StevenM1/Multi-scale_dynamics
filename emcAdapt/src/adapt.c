#include <math.h>
#include <Rmath.h>
#include <stdlib.h>
#include <R.h>
#include <stdbool.h>
#include <stdio.h>

// delta rule
void adaptDelta(int *nTrials,
                int *nChoices,
                double *values,
                double *adaptedValues, // of length ntrials*n_pairs
                double *predictionErrors,
                double *outcome,
                double *learningRates) {
  // nTrials: number of trials
  // nChoices: total number of choice options. For example, 3 pairs of 2 stimuli = 6 choices
  // values: initial value of each choice option (length nChoices)
  // VV: output for trial-by-trial values of each choice option. In R, a matrix of size (nTrials, nChoices); here in C, an array of length (nTrials*nChoices)
  // PE: output for trial-by-trial prediction errors for each choie option. Identical size as VV
  // outcome: (nTrials, nChoices): trial by trial outcomes per choice option. NA is no output.
  // eta1, eta2: floats with the learning rates (positive and negative, respectively)

  // declare some variables
  //  double updateValue;
  double this_lr;
  int mat_idx = 0;
  int nt = *nTrials;
  int nc = *nChoices;
  static double dv[64] = {0};   // WARNING: MAXIMUM NUMBER OF CHOICES IS FIXED HERE TO 64. Anything more will crash.

  // Loop over trials i
  for(unsigned int i = 0; i < nt; i++) {
    //    printf("Trial N: %d\n", i);
    //    this_eta1 = eta1[i];  // trialwise learning rate for positive PE
    //    this_eta2 = eta2[i];  // trialwise learning rate for negative PE

    // Loop over choice alternatives
    for(unsigned int ch = 0; ch < nc; ch++) {
      //      printf("Choice option: %d\n", ch);
      mat_idx = nt*ch+i;  // Where in the VV-matrix and outcome-matrix are we?
      adaptedValues[mat_idx] = values[ch]; // Offload current values

      // If the outcome for this choice was not NA, update (note that o_ == o_ returns FALSE if o_ is NA)
      if(outcome[mat_idx] == outcome[mat_idx]) {
        this_lr = learningRates[mat_idx];
        //        printf("Outcome is: %.3f, this is NOT NA!\n", outcome[mat_idx]);
        dv[ch] = outcome[mat_idx] - values[ch];  // prediction error dv = outcome choice - value of choice
        predictionErrors[mat_idx] = dv[ch];  // offload PE

        // Do we update with eta1 or eta2?
        //        updateValue = dv[ch] > 0 ? this_eta1 : this_eta2;  // Ternary expression  (if-else in a one-liner)

        values[ch] = values[ch] + this_lr*dv[ch];  // Update value
      }
    }
  }
}

// // special updating for thresholds based on drift
// void adaptThresholds(int *nTrials,
//                      double *b0,
//                      double *weight,
//                      //                     int *nChoices,
//                      double *values,
//                      double *adaptedValues,
//                      double *predictionErrors,
//                      double *rts,
//                      double *learningRates) {
//   // nTrials: number of trials
//   // declare some variables
//   double this_lr;
//   int mat_idx = 0;
//   int nt = *nTrials;
//   //  double b0_ = *b0;
//   double weight_ = *weight;
//   double dv;
//
//   // Loop over trials i
//   for(unsigned int i = 0; i < nt; i++) {
//     //    printf("Trial N: %d\n", i);
//     //    this_eta1 = eta1[i];  // trialwise learning rate for positive PE
//     //    this_eta2 = eta2[i];  // trialwise learning rate for negative PE
//
//     // Loop over choice alternatives
//     //    for(unsigned int ch = 0; ch < nc; ch++) {
//     //      printf("Choice option: %d\n", ch);
//     //      mat_idx = nt*ch+i;  // Where in the VV-matrix and outcome-matrix are we?
//     for(unsigned int ch = 0; ch < 2; ch++) {
//       mat_idx = nt*ch+i;
//       adaptedValues[mat_idx] = values[ch]; // Offload current values
//     }
//     this_lr = learningRates[i];
//     // prediction error here is defined as (b/rt - q); values[0] = b, values[1] = q
//     dv = (values[0]/rts[i]) - values[1];
//     //printf("values: %f\n", values[0]);
//
//
//     //    values[2] <- values[2] + learningRates[trial,2]*dv
//     values[1] = values[1] + this_lr*dv;   // update q
//     values[0] = b0[i] + weight_*values[1];  // update threshold
//
//     predictionErrors[mat_idx] = dv;
//   }
// }

// special updating for thresholds based on drift
void adaptThresholds(int *nTrials,
                     double *b0,
                     double *weight,
                     double *values,
                     double *adaptedValues,
                     double *predictionErrors,
                     double *rts,
                     double *learningRates) {
  // nTrials: number of trials
  // declare some variables
  double this_lr;
  int mat_idx = 0;
  int nt = *nTrials;
  double weight_ = *weight;
  double dv;

  // Loop over trials i
  for(unsigned int i = 0; i < nt; i++) {
    // if trial > 0, save last trial's data
    if(i > 0) {
      // for(unsigned int ch = 0; ch < 2; ch++) {
      //   mat_idx = nt*ch+(i-1);  // PREVIOUS TRIAL
      //   adaptedValues[mat_idx] = values[ch]; // Offload current values
      // }
      //  // and update Q-values
      //  values[1] = values[1] + this_lr*dv;
    }
    this_lr = learningRates[i];
    // prediction error here is defined as (b/rt - q); values[0] = b, values[1] = q
    values[0] = b0[i] + weight_*values[1];  // update threshold
    dv = (values[0]/rts[i]) - values[1];    // calculate prediction error
    predictionErrors[i] = dv;               // save prediction error

    // update channels
    for(unsigned int ch = 0; ch < 2; ch++) {
      mat_idx = nt*ch+i;  // find indx in adaptedValues
      adaptedValues[mat_idx] = values[ch]; // Offload current values
    }
    // and update Q-values
    values[1] = values[1] + this_lr*dv;
  }
}


// delta rule
void adaptDelta2LR(int *nTrials,
                   int *nChoices,
                   double *values,
                   double *adaptedValues, // of length ntrials*n_pairs
                   double *predictionErrors,
                   double *outcome,
                   double *learningRatesPos,
                   double *learningRatesNeg) {
  // nTrials: number of trials
  // nChoices: total number of choice options. For example, 3 pairs of 2 stimuli = 6 choices
  // values: initial value of each choice option (length nChoices)
  // VV: output for trial-by-trial values of each choice option. In R, a matrix of size (nTrials, nChoices); here in C, an array of length (nTrials*nChoices)
  // PE: output for trial-by-trial prediction errors for each choie option. Identical size as VV
  // outcome: (nTrials, nChoices): trial by trial outcomes per choice option. NA is no output.
  // learningRatesPos, learningRatesNeg: floats with the learning rates (positive and negative, respectively)

  // declare some variables
  //  double updateValue;
  double this_lr;
  int mat_idx = 0;
  int nt = *nTrials;
  int nc = *nChoices;
  static double dv[64] = {0};   // WARNING: MAXIMUM NUMBER OF CHOICES IS FIXED HERE TO 64. Anything more will crash.

  // Loop over trials i
  for(unsigned int i = 0; i < nt; i++) {
    //    printf("Trial N: %d\n", i);
    //    this_eta1 = eta1[i];  // trialwise learning rate for positive PE
    //    this_eta2 = eta2[i];  // trialwise learning rate for negative PE

    // Loop over choice alternatives
    for(unsigned int ch = 0; ch < nc; ch++) {
      //      printf("Choice option: %d\n", ch);
      mat_idx = nt*ch+i;  // Where in the VV-matrix and outcome-matrix are we?
      adaptedValues[mat_idx] = values[ch]; // Offload current values

      // If the outcome for this choice was not NA, update (note that o_ == o_ returns FALSE if o_ is NA)
      if(outcome[mat_idx] == outcome[mat_idx]) {
//        this_lr = learningRates[mat_idx];
        //        printf("Outcome is: %.3f, this is NOT NA!\n", outcome[mat_idx]);
        dv[ch] = outcome[mat_idx] - values[ch];  // prediction error dv = outcome choice - value of choice
        predictionErrors[mat_idx] = dv[ch];  // offload PE

        // Do we update with a positive or a negative learning rate?
        this_lr = dv[ch] >= 0 ? learningRatesPos[mat_idx] : learningRatesNeg[mat_idx];  // Ternary expression  (if-else in a one-liner)

        values[ch] = values[ch] + this_lr*dv[ch];  // Update value
      }
    }
  }
}


// Volatile Kalman filer (linear)
void adaptVKF(int *nTrials,
              int *nChoices,

              double *predictions,               // m (size nc)
              double *adaptedPredictions,        // m (size nt x nc)
              double *predictionErrors,          // delta_m (size nt x nc)
              double *learningRates,             // k (size nt x nc)

              double *volatilities,               // v (size nc)
              double *adaptedVolatilities,        // v (size nt x nc)
              double *volatilityPredictionErrors, // delta_v (size nt x nc)
              double *volatilityLearningRates,    // lambda (size nt x nc)

              double *uncertainties,              // w (size nc)
              double *adaptedUncertainties,       // w (size nt x nc)
              // double *sigma2,                     // size nc (never updated)

              double *outcomes                   // rewards size (nt x nc)
              ) {

  // nTrials: number of trials
  // nChoices: total number of choice options. For example, 3 pairs of 2 stimuli = 6 choices
  // predictions: initial value of each choice option (length nChoices)
  // adaptedPredictions: output for trial-by-trial values of each choice option. In R, a matrix of size (nTrials, nChoices); here in C, an array of length (nTrials*nChoices)
  // predictionErrors: output for trial-by-trial prediction errors for each chocie option. Identical size as adaptedPredictions
  // learningRates: output for trial-by-trial learning rates per choice option

  // volatilities: initial value of each choice option's volatility
  // adaptedVolatility: output for trial-by-trial volatility of each choice option
  // volatilityPredictionErrors: output for trial-by-trial prediction errors for each choice option's volatility
  // volatilityLearningRates: input for trial-by-trial volatility learning rates per choice option

  // uncertainties: initial value of each choice option's uncertainty
  // adaptedUncertainties: output for uncertainties

  // outcomes: (nTrials, nChoices): trial by trial outcomes per choice option. NA is no output.

  // declare some variables

  // placeholders
  double this_lr;             // k
  double this_volatility;     // v
  double this_uncertainty;    // w
  double this_volatility_lr;  // lambda

  double mpre;
  double wpre;
  double wcov;
  double delta_v;

  // indexers, sizes
  int mat_idx = 0;
  int nt = *nTrials;
  int nc = *nChoices;

  // to keep track of prediction errors per choice option
  static double delta_prediction[64] = {0};              // WARNING: MAXIMUM NUMBER OF CHOICES IS FIXED HERE TO 64. Anything more will crash.
  static double delta_volatility[64] = {0};   // WARNING: MAXIMUM NUMBER OF CHOICES IS FIXED HERE TO 64. Anything more will crash.

  // copy initial uncertainties to sigma2
  static double sigma2[64] = {0};  // sigma2 = initial uncertainty, and static
  for(unsigned int ch = 0; ch < nc; ch++) {
    sigma2[ch] = uncertainties[ch];
  }

  // Loop over trials i
  for(unsigned int i = 0; i < nt; i++) {
    //    printf("Trial N: %d\n", i);
    //    this_eta1 = eta1[i];  // trialwise learning rate for positive PE
    //    this_eta2 = eta2[i];  // trialwise learning rate for negative PE

    // Loop over choice alternatives
    for(unsigned int ch = 0; ch < nc; ch++) {
      //      printf("Choice option: %d\n", ch);
      mat_idx = nt*ch+i;  // Where in the VV-matrix and outcome-matrix are we?
      adaptedPredictions[mat_idx] = predictions[ch];     // Offload current values
      adaptedVolatilities[mat_idx] = volatilities[ch];   // Offload current values
      adaptedUncertainties[mat_idx] = uncertainties[ch]; // Offload current values

      // If the outcome for this choice was not NA, update (note that o_ == o_ returns FALSE if o_ is NA)
      if(outcomes[mat_idx] == outcomes[mat_idx]) {

        this_volatility_lr = volatilityLearningRates[mat_idx];
        mpre = predictions[ch];
        wpre = uncertainties[ch];

        // printf("Sigma2: %f", sigma2[ch]);
        // printf("  , Uncertainty: %f\n", uncertainties[ch]);

        // prediction error dv = outcome choice - value of choice,
        delta_prediction[ch] = outcomes[mat_idx] - predictions[ch];

        // learning rate,    Eq 9
        this_lr = (uncertainties[ch] + volatilities[ch]) / (uncertainties[ch] + volatilities[ch] + sigma2[ch]);
        // printf("uncertainties[ch]: %f, volatilities[ch]: %f, sigma2[ch]: %f, resulting this_lr: %f\n", uncertainties[ch], volatilities[ch], sigma2[ch], this_lr);

        // prediction, eq 10
        predictions[ch] = predictions[ch] + this_lr * delta_prediction[ch];
        // uncertainty, Eq 11
        uncertainties[ch] = (1 - this_lr) * (uncertainties[ch] + volatilities[ch]);

        // covariance of uncertainty? Eq 12
        wcov = (1 - this_lr) * wpre;
        // printf("this_lr: %f, wpre: %f, resultig wcov: %f\n", this_lr, wpre, wcov);

        // volatility error // FROM https://github.com/payampiray/VKF/blob/master/vkf.m !!!
        delta_volatility[ch] =  pow(predictions[ch] - mpre, 2) + uncertainties[ch] + wpre - 2 * wcov - volatilities[ch];
        // printf("predictions[ch]: %f, mpre: %f, uncertainties[ch]: %f, wpre: %f, wcov: %f, volatilities[ch]: %f, resulting delta_volatility[ch]: %f\n", predictions[ch], mpre, uncertainties[ch], wpre, wcov, volatilities[ch], delta_volatility[ch]);

        // update volatility, Eq 13
        volatilities[ch] = volatilities[ch] + this_volatility_lr*delta_volatility[ch];
        // printf("volatilities[ch]: %f = volatilities[ch]: %f + this_volatility_lr[ch]: %f * delta_volatility: %f\n", volatilities[ch], volatilities[ch], this_volatility_lr, delta_volatility[ch]);

        // printf("\n");
        // save what we want to save: PEs, volatility PEs, learning rates
        predictionErrors[mat_idx] = delta_prediction[ch];  // offload PE
        learningRates[mat_idx] = this_lr;
        volatilityPredictionErrors[mat_idx] = delta_volatility[ch];
      }
    }
  }
}


// Volatile Kalman filer (binary)
void adaptVKFbinary(int *nTrials,
              int *nChoices,

              double *predictions,               // m (size nc)
              double *adaptedPredictions,        // m (size nt x nc)
              double *predictionErrors,          // delta_m (size nt x nc)
              double *learningRates,             // k (size nt x nc)

              double *volatilities,               // v (size nc)
              double *adaptedVolatilities,        // v (size nt x nc)
              double *volatilityPredictionErrors, // delta_v (size nt x nc)
              double *volatilityLearningRates,    // lambda (size nt x nc)

              double *uncertainties,              // w (size nc)
              double *adaptedUncertainties,       // w (size nt x nc)
              // double *sigma2,                     // size nc (never updated)

              double *outcomes                   // rewards size (nt x nc)
) {

  // nTrials: number of trials
  // nChoices: total number of choice options. For example, 3 pairs of 2 stimuli = 6 choices
  // predictions: initial value of each choice option (length nChoices)
  // adaptedPredictions: output for trial-by-trial values of each choice option. In R, a matrix of size (nTrials, nChoices); here in C, an array of length (nTrials*nChoices)
  // predictionErrors: output for trial-by-trial prediction errors for each chocie option. Identical size as adaptedPredictions
  // learningRates: output for trial-by-trial learning rates per choice option

  // volatilities: initial value of each choice option's volatility
  // adaptedVolatility: output for trial-by-trial volatility of each choice option
  // volatilityPredictionErrors: output for trial-by-trial prediction errors for each choice option's volatility
  // volatilityLearningRates: input for trial-by-trial volatility learning rates per choice option

  // uncertainties: initial value of each choice option's uncertainty
  // adaptedUncertainties: output for uncertainties

  // outcomes: (nTrials, nChoices): trial by trial outcomes per choice option. NA is no output.

  // declare some variables

  // placeholders
  double this_lr;             // alpha
  double this_kalman_gain;    // k
  double this_volatility;     // v
  double this_uncertainty;    // w
  double this_volatility_lr;  // lambda

  double mpre;
  double wpre;
  double wcov;
  double delta_v;

  // indexers, sizes
  int mat_idx = 0;
  int nt = *nTrials;
  int nc = *nChoices;

  // to keep track of prediction errors per choice option
  static double delta_prediction[64] = {0};              // WARNING: MAXIMUM NUMBER OF CHOICES IS FIXED HERE TO 64. Anything more will crash.
  static double delta_volatility[64] = {0};   // WARNING: MAXIMUM NUMBER OF CHOICES IS FIXED HERE TO 64. Anything more will crash.

  // copy initial uncertainties to omega
  static double omega[64] = {0};  // omega = initial uncertainty, and static
  for(unsigned int ch = 0; ch < nc; ch++) {
    omega[ch] = uncertainties[ch];
  }

  // Loop over trials i
  for(unsigned int i = 0; i < nt; i++) {
    //    printf("Trial N: %d\n", i);
    //    this_eta1 = eta1[i];  // trialwise learning rate for positive PE
    //    this_eta2 = eta2[i];  // trialwise learning rate for negative PE

    // Loop over choice alternatives
    for(unsigned int ch = 0; ch < nc; ch++) {
      //      printf("Choice option: %d\n", ch);
      mat_idx = nt*ch+i;  // Where in the VV-matrix and outcome-matrix are we?
      adaptedPredictions[mat_idx] = predictions[ch];     // Offload current values
      adaptedVolatilities[mat_idx] = volatilities[ch];   // Offload current values
      adaptedUncertainties[mat_idx] = uncertainties[ch]; // Offload current values

      // If the outcome for this choice was not NA, update (note that o_ == o_ returns FALSE if o_ is NA)
      if(outcomes[mat_idx] == outcomes[mat_idx]) {

        this_volatility_lr = volatilityLearningRates[mat_idx];
        mpre = predictions[ch];
        wpre = uncertainties[ch];

        // printf("Uncertainty: %f\n", sigma2[ch]);

        // prediction error dv = outcome choice - value of choice, eq 16.1
        delta_prediction[ch] = outcomes[mat_idx] - (1 / (1+exp(-predictions[ch])));

        // kalman gain,    Eq 14
        this_kalman_gain = (uncertainties[ch] + volatilities[ch]) / (uncertainties[ch] + volatilities[ch] + omega[ch]);
        // printf("uncertainties[ch]: %f, volatilities[ch]: %f, sigma2[ch]: %f, resulting this_lr: %f\n", uncertainties[ch], volatilities[ch], sigma2[ch], this_lr);

        // learning rate, Eq 15
        this_lr = sqrt(uncertainties[ch] + volatilities[ch]);

        // prediction, eq 16.2
        predictions[ch] = predictions[ch] + this_lr * delta_prediction[ch];
        // uncertainty, Eq 11
        uncertainties[ch] = (1 - this_kalman_gain) * (uncertainties[ch] + volatilities[ch]);

        // covariance of uncertainty? Eq 12
        wcov = (1 - this_kalman_gain) * wpre;
        // printf("this_kalman_gain: %f, wpre: %f, resultig wcov: %f\n", this_kalman_gain, wpre, wcov);

        // volatility error FROM https://github.com/payampiray/VKF/blob/master/vkf_bin.m !!!
        delta_volatility[ch] =  pow(predictions[ch] - mpre, 2) + uncertainties[ch] + wpre - 2 * wcov - volatilities[ch];
        // printf("predictions[ch]: %f, mpre: %f, uncertainties[ch]: %f, wpre: %f, wcov: %f, volatilities[ch]: %f, resulting delta_volatility[ch]: %f\n", predictions[ch], mpre, uncertainties[ch], wpre, wcov, volatilities[ch], delta_volatility[ch]);

        // update volatility, Eq 13
        volatilities[ch] = volatilities[ch] + this_volatility_lr*delta_volatility[ch];
        // printf("volatilities[ch]: %f = volatilities[ch]: %f + this_volatility_lr[ch]: %f * delta_volatility: %f\n", volatilities[ch], volatilities[ch], this_volatility_lr, delta_volatility[ch]);

        // printf("\n");
        // save what we want to save: PEs, volatility PEs, learning rates
        predictionErrors[mat_idx] = delta_prediction[ch];  // offload PE
        learningRates[mat_idx] = this_lr;
        volatilityPredictionErrors[mat_idx] = delta_volatility[ch];
      }
    }
  }
}

/// different implementations:
// Volatile Kalman filer (linear)
void adaptVKF2(int *nTrials,
              int *nChoices,

              double *predictions,               // m (size nc)
              double *adaptedPredictions,        // m (size nt x nc)
              double *predictionErrors,          // delta_m (size nt x nc)
              double *learningRates,             // k (size nt x nc)

              double *volatilities,               // v (size nc)
              double *adaptedVolatilities,        // v (size nt x nc)
              double *volatilityPredictionErrors, // delta_v (size nt x nc)
              double *volatilityLearningRates,    // lambda (size nt x nc)

              double *uncertainties,              // w (size nc)
              double *adaptedUncertainties,       // w (size nt x nc)
              // double *sigma2,                     // size nc (never updated)

              double *outcomes                   // rewards size (nt x nc)
) {

  // nTrials: number of trials
  // nChoices: total number of choice options. For example, 3 pairs of 2 stimuli = 6 choices
  // predictions: initial value of each choice option (length nChoices)
  // adaptedPredictions: output for trial-by-trial values of each choice option. In R, a matrix of size (nTrials, nChoices); here in C, an array of length (nTrials*nChoices)
  // predictionErrors: output for trial-by-trial prediction errors for each chocie option. Identical size as adaptedPredictions
  // learningRates: output for trial-by-trial learning rates per choice option

  // volatilities: initial value of each choice option's volatility
  // adaptedVolatility: output for trial-by-trial volatility of each choice option
  // volatilityPredictionErrors: output for trial-by-trial prediction errors for each choice option's volatility
  // volatilityLearningRates: input for trial-by-trial volatility learning rates per choice option

  // uncertainties: initial value of each choice option's uncertainty
  // adaptedUncertainties: output for uncertainties

  // outcomes: (nTrials, nChoices): trial by trial outcomes per choice option. NA is no output.

  // declare some variables

  // placeholders
  double this_lr;             // k
  double this_volatility;     // v
  double this_uncertainty;    // w
  double this_volatility_lr;  // lambda

  double mpre;
  double wpre;
  double wcov;
  double delta_v;

  // indexers, sizes
  int mat_idx = 0;
  int nt = *nTrials;
  int nc = *nChoices;

  // to keep track of prediction errors per choice option
  static double delta_prediction[64] = {0};              // WARNING: MAXIMUM NUMBER OF CHOICES IS FIXED HERE TO 64. Anything more will crash.
  static double delta_volatility[64] = {0};   // WARNING: MAXIMUM NUMBER OF CHOICES IS FIXED HERE TO 64. Anything more will crash.

  // copy initial uncertainties to sigma2
  static double sigma2[64] = {0};  // sigma2 = initial uncertainty, and static
  for(unsigned int ch = 0; ch < nc; ch++) {
    sigma2[ch] = uncertainties[ch];
  }

  // Loop over trials i
  for(unsigned int i = 0; i < nt; i++) {
    //    printf("Trial N: %d\n", i);
    //    this_eta1 = eta1[i];  // trialwise learning rate for positive PE
    //    this_eta2 = eta2[i];  // trialwise learning rate for negative PE

    // Loop over choice alternatives
    for(unsigned int ch = 0; ch < nc; ch++) {
      //      printf("Choice option: %d\n", ch);
      mat_idx = nt*ch+i;  // Where in the VV-matrix and outcome-matrix are we?
      adaptedPredictions[mat_idx] = predictions[ch];     // Offload current values
      adaptedVolatilities[mat_idx] = volatilities[ch];   // Offload current values
      adaptedUncertainties[mat_idx] = uncertainties[ch]; // Offload current values

      // If the outcome for this choice was not NA, update (note that o_ == o_ returns FALSE if o_ is NA)
      if(outcomes[mat_idx] == outcomes[mat_idx]) {

        this_volatility_lr = volatilityLearningRates[mat_idx];
        mpre = predictions[ch];
        wpre = uncertainties[ch];

        // printf("Sigma2: %f", sigma2[ch]);
        // printf("  , Uncertainty: %f\n", uncertainties[ch]);

        // prediction error dv = outcome choice - value of choice,
        delta_prediction[ch] = outcomes[mat_idx] - predictions[ch];

        // learning rate,    Eq 9
        this_lr = (uncertainties[ch] + volatilities[ch]) / (uncertainties[ch] + volatilities[ch] + sigma2[ch]);
        // printf("uncertainties[ch]: %f, volatilities[ch]: %f, sigma2[ch]: %f, resulting this_lr: %f\n", uncertainties[ch], volatilities[ch], sigma2[ch], this_lr);

        // prediction, eq 10
        predictions[ch] = predictions[ch] + this_lr * delta_prediction[ch];
        // uncertainty, Eq 11
        uncertainties[ch] = (1 - this_lr) * (uncertainties[ch] + volatilities[ch]);

        // covariance of uncertainty? Eq 12
//        wcov = (1 - this_lr) * wpre;
        // printf("this_lr: %f, wpre: %f, resultig wcov: %f\n", this_lr, wpre, wcov);

        // volatility error  FROM https://github.com/payampiray/piray_daw_2020_ploscb/blob/master/vkf_lin.m !!!
        delta_volatility[ch] = pow(this_lr, 2) * pow(delta_prediction[ch], 2) + this_lr*wpre - this_lr*volatilities[ch];
       // delta_volatility[ch] =  pow(predictions[ch] - mpre, 2) + uncertainties[ch] + wpre - 2 * wcov - volatilities[ch];
        // printf("predictions[ch]: %f, mpre: %f, uncertainties[ch]: %f, wpre: %f, wcov: %f, volatilities[ch]: %f, resulting delta_volatility[ch]: %f\n", predictions[ch], mpre, uncertainties[ch], wpre, wcov, volatilities[ch], delta_volatility[ch]);

        // update volatility, Eq 13
        volatilities[ch] = volatilities[ch] + this_volatility_lr*delta_volatility[ch];
        // printf("volatilities[ch]: %f = volatilities[ch]: %f + this_volatility_lr[ch]: %f * delta_volatility: %f\n", volatilities[ch], volatilities[ch], this_volatility_lr, delta_volatility[ch]);

        // printf("\n");
        // save what we want to save: PEs, volatility PEs, learning rates
        predictionErrors[mat_idx] = delta_prediction[ch];  // offload PE
        learningRates[mat_idx] = this_lr;
        volatilityPredictionErrors[mat_idx] = delta_volatility[ch];
      }
    }
  }
}


// Volatile Kalman filer (binary)
void adaptVKFbinary2(int *nTrials,
                    int *nChoices,

                    double *predictions,               // m (size nc)
                    double *adaptedPredictions,        // m (size nt x nc)
                    double *predictionErrors,          // delta_m (size nt x nc)
                    double *learningRates,             // k (size nt x nc)

                    double *volatilities,               // v (size nc)
                    double *adaptedVolatilities,        // v (size nt x nc)
                    double *volatilityPredictionErrors, // delta_v (size nt x nc)
                    double *volatilityLearningRates,    // lambda (size nt x nc)

                    double *uncertainties,              // w (size nc)
                    double *adaptedUncertainties,       // w (size nt x nc)
                    // double *sigma2,                     // size nc (never updated)

                    double *outcomes                   // rewards size (nt x nc)
) {

  // nTrials: number of trials
  // nChoices: total number of choice options. For example, 3 pairs of 2 stimuli = 6 choices
  // predictions: initial value of each choice option (length nChoices)
  // adaptedPredictions: output for trial-by-trial values of each choice option. In R, a matrix of size (nTrials, nChoices); here in C, an array of length (nTrials*nChoices)
  // predictionErrors: output for trial-by-trial prediction errors for each chocie option. Identical size as adaptedPredictions
  // learningRates: output for trial-by-trial learning rates per choice option

  // volatilities: initial value of each choice option's volatility
  // adaptedVolatility: output for trial-by-trial volatility of each choice option
  // volatilityPredictionErrors: output for trial-by-trial prediction errors for each choice option's volatility
  // volatilityLearningRates: input for trial-by-trial volatility learning rates per choice option

  // uncertainties: initial value of each choice option's uncertainty
  // adaptedUncertainties: output for uncertainties

  // outcomes: (nTrials, nChoices): trial by trial outcomes per choice option. NA is no output.

  // declare some variables

  // placeholders
  double this_lr;             // alpha
  double this_kalman_gain;    // k
  double this_volatility;     // v
  double this_uncertainty;    // w
  double this_volatility_lr;  // lambda

  double mpre;
  double wpre;
  double wcov;
  double delta_v;

  // indexers, sizes
  int mat_idx = 0;
  int nt = *nTrials;
  int nc = *nChoices;

  // to keep track of prediction errors per choice option
  static double delta_prediction[64] = {0};              // WARNING: MAXIMUM NUMBER OF CHOICES IS FIXED HERE TO 64. Anything more will crash.
  static double delta_volatility[64] = {0};   // WARNING: MAXIMUM NUMBER OF CHOICES IS FIXED HERE TO 64. Anything more will crash.

  // copy initial uncertainties to omega
  static double omega[64] = {0};  // omega = initial uncertainty, and static
  for(unsigned int ch = 0; ch < nc; ch++) {
    omega[ch] = uncertainties[ch];
  }

  // Loop over trials i
  for(unsigned int i = 0; i < nt; i++) {
    //    printf("Trial N: %d\n", i);
    //    this_eta1 = eta1[i];  // trialwise learning rate for positive PE
    //    this_eta2 = eta2[i];  // trialwise learning rate for negative PE

    // Loop over choice alternatives
    for(unsigned int ch = 0; ch < nc; ch++) {
      //      printf("Choice option: %d\n", ch);
      mat_idx = nt*ch+i;  // Where in the VV-matrix and outcome-matrix are we?
      adaptedPredictions[mat_idx] = predictions[ch];     // Offload current values
      adaptedVolatilities[mat_idx] = volatilities[ch];   // Offload current values
      adaptedUncertainties[mat_idx] = uncertainties[ch]; // Offload current values

      // If the outcome for this choice was not NA, update (note that o_ == o_ returns FALSE if o_ is NA)
      if(outcomes[mat_idx] == outcomes[mat_idx]) {

        this_volatility_lr = volatilityLearningRates[mat_idx];
        mpre = predictions[ch];
        wpre = uncertainties[ch];

        // printf("Uncertainty: %f\n", sigma2[ch]);

        // prediction error dv = outcome choice - value of choice, eq 16.1
        delta_prediction[ch] = outcomes[mat_idx] - (1 / (1+exp(-predictions[ch])));

        // kalman gain,    Eq 14
        this_kalman_gain = (uncertainties[ch] + volatilities[ch]) / (uncertainties[ch] + volatilities[ch] + omega[ch]);
        // printf("uncertainties[ch]: %f, volatilities[ch]: %f, sigma2[ch]: %f, resulting this_lr: %f\n", uncertainties[ch], volatilities[ch], sigma2[ch], this_lr);

        // learning rate, Eq 15
        this_lr = sqrt(uncertainties[ch] + volatilities[ch]);

        // prediction, eq 16.2
        predictions[ch] = predictions[ch] + this_lr * delta_prediction[ch];
        // uncertainty, Eq 11
        uncertainties[ch] = (1 - this_kalman_gain) * (uncertainties[ch] + volatilities[ch]);

        // covariance of uncertainty? Eq 12
//        wcov = (1 - this_kalman_gain) * wpre;
        // printf("this_kalman_gain: %f, wpre: %f, resultig wcov: %f\n", this_kalman_gain, wpre, wcov);

        // volatility error FROM https://github.com/payampiray/piray_daw_2020_ploscb/blob/master/vkf_bin.m !!!
        //v           = v +lambda.*(delta.^2 + k.*wpre - k.*v);
        delta_volatility[ch] = pow(this_lr*delta_prediction[ch], 2) + this_kalman_gain*wpre - this_kalman_gain*volatilities[ch];
//        delta_volatility[ch] =  pow(predictions[ch] - mpre, 2) + uncertainties[ch] + wpre - 2 * wcov - volatilities[ch];
//        printf("predictions[ch]: %f, mpre: %f, uncertainties[ch]: %f, wpre: %f, wcov: %f, volatilities[ch]: %f, resulting delta_volatility[ch]: %f\n", predictions[ch], mpre, uncertainties[ch], wpre, wcov, volatilities[ch], delta_volatility[ch]);

        // update volatility, Eq 13
        volatilities[ch] = volatilities[ch] + this_volatility_lr*delta_volatility[ch];
        // printf("volatilities[ch]: %f = volatilities[ch]: %f + this_volatility_lr[ch]: %f * delta_volatility: %f\n", volatilities[ch], volatilities[ch], this_volatility_lr, delta_volatility[ch]);

        // printf("\n");
        // save what we want to save: PEs, volatility PEs, learning rates
        predictionErrors[mat_idx] = delta_prediction[ch];  // offload PE
        learningRates[mat_idx] = this_lr;
        volatilityPredictionErrors[mat_idx] = delta_volatility[ch];
      }
    }
  }
}




// // Dual learning rates
// void adaptSARSA2LR(int *nTrials,
//                    int *nChoices,
//                    double *values,
//                    double *adaptedValues, // of length ntrials*n_pairs
//                    double *predictionErrors,
//                    double *outcome,
//                    double *learningRatesPos,
//                    double *learningRatesNeg) {
//   // nTrials: number of trials
//   // nChoices: total number of choice options. For example, 3 pairs of 2 stimuli = 6 choices
//   // values: initial value of each choice option (length nChoices)
//   // VV: output for trial-by-trial values of each choice option. In R, a matrix of size (nTrials, nChoices); here in C, an array of length (nTrials*nChoices)
//   // PE: output for trial-by-trial prediction errors for each choie option. Identical size as VV
//   // outcome: (nTrials, nChoices): trial by trial outcomes per choice option. NA is no output.
//   // eta1, eta2: floats with the learning rates (positive and negative, respectively)
//
//   // declare some variables
//   double updateValue;
//   double this_lr;
//   int mat_idx = 0;
//   int nt = *nTrials;
//   int nc = *nChoices;
//   static double dv[64] = {0};   // WARNING: MAXIMUM NUMBER OF CHOICES IS FIXED HERE TO 64. Anything more will crash.
//
//   // Loop over trials i
//   for(unsigned int i = 0; i < nt; i++) {
//     //    printf("Trial N: %d\n", i);
//     //    this_eta1 = eta1[i];  // trialwise learning rate for positive PE
//     //    this_eta2 = eta2[i];  // trialwise learning rate for negative PE
//
//     // Loop over choice alternatives
//     for(unsigned int ch = 0; ch < nc; ch++) {
//       //      printf("Choice option: %d\n", ch);
//       mat_idx = nt*ch+i;  // Where in the VV-matrix and outcome-matrix are we?
//       adaptedValues[mat_idx] = values[ch]; // Offload current values
//
//       // If the outcome for this choice was not NA, update (note that o_ == o_ returns FALSE if o_ is NA)
//       if(outcome[mat_idx] == outcome[mat_idx]) {
//         //        this_lr = learningRates[mat_idx];
//         //        printf("Outcome is: %.3f, this is NOT NA!\n", outcome[mat_idx]);
//         dv[ch] = outcome[mat_idx] - values[ch];  // prediction error dv = outcome choice - value of choice
//         predictionErrors[mat_idx] = dv[ch];  // offload PE
//
//         // Do we update with eta1 or eta2?
//         this_lr = dv[ch] > 0 ? learningRatesPos[mat_idx] : learningRatesNeg[mat_idx];  // Ternary expression  (if-else in a one-liner)
//
//         values[ch] = values[ch] + this_lr*dv[ch];  // Update value
//       }
//     }
//   }
// }
//
//
//
// void adaptSARSARisk(int *nTrials,
//                      int *nChoices,
//                      double *values,
//                      double *adaptedValues, // of length ntrials*n_pairs
//                      double *riskValues,
//                      double *adaptedRiskValues,
//                      double *predictionErrors,
//                      double *riskPredictionErrors,
//                      double *outcome,
//                      double *learningRates,
//                      double *riskLearningRates) {
//   // nTrials: number of trials
//   // nChoices: total number of choice options. For example, 3 pairs of 2 stimuli = 6 choices
//   // values: initial value of each choice option (length nChoices)
//   // VV: output for trial-by-trial values of each choice option. In R, a matrix of size (nTrials, nChoices); here in C, an array of length (nTrials*nChoices)
//   // PE: output for trial-by-trial prediction errors for each choie option. Identical size as VV
//   // outcome: (nTrials, nChoices): trial by trial outcomes per choice option. NA is no output.
//   // eta1, eta2: floats with the learning rates (positive and negative, respectively)
//
//   // declare some variables
//   //  double updateValue;
//   double this_lr;
//   double this_risk_lr;
//   int mat_idx = 0;
//   int nt = *nTrials;
//   int nc = *nChoices;
//   static double dv[64] = {0};   // WARNING: MAXIMUM NUMBER OF CHOICES IS FIXED HERE TO 64. Anything more will crash.
//   static double riskdv[64] = {0};
//
//   // Loop over trials i
//   for(unsigned int i = 0; i < nt; i++) {
//     //    printf("Trial N: %d\n", i);
//     //    this_eta1 = eta1[i];  // trialwise learning rate for positive PE
//     //    this_eta2 = eta2[i];  // trialwise learning rate for negative PE
//
//     // Loop over choice alternatives
//     for(unsigned int ch = 0; ch < nc; ch++) {
//       //      printf("Choice option: %d\n", ch);
//       mat_idx = nt*ch+i;  // Where in the VV-matrix and outcome-matrix are we?
//       adaptedValues[mat_idx] = values[ch]; // Offload current values
//       adaptedRiskValues[mat_idx] = riskValues[ch]; // offload risk
//
//       // If the outcome for this choice was not NA, update (note that o_ == o_ returns FALSE if o_ is NA)
//       if(outcome[mat_idx] == outcome[mat_idx]) {
//         this_lr = learningRates[mat_idx];
//         this_risk_lr = learningRates[mat_idx];
//
//         // Compute reward prediction error
//         dv[ch] = outcome[mat_idx] - values[ch];  // prediction error dv = outcome choice - value of choice
//         predictionErrors[mat_idx] = dv[ch];  // offload PE
//
//         // Compute risk prediction error
//         riskdv[ch] = dv[ch]*dv[ch] - riskValues[ch];
//         riskPredictionErrors[mat_idx] = riskdv[ch];
//
//         // Update risk
//         riskValues[ch] = riskValues[ch] + this_risk_lr * riskdv[ch];
//
//         // Update reward expectation, scale PE by risk
//         values[ch] = values[ch] + this_lr*(dv[ch] / sqrt(riskValues[ch]));  // Update value
//
//         // offload actual learning rate
// //        learningRates[mat_idx] = this_variance*this_lr;
//       }
//     }
//   }
// }
//
//
// void adaptSARSAvarRL(int *nTrials,
//                     int *nChoices,
//                     double *values,
//                     double *adaptedValues, // of length ntrials*n_pairs
//                     double *predictionErrors,
//                     double *outcome,
//                     double *learningRates) {
//   // nTrials: number of trials
//   // nChoices: total number of choice options. For example, 3 pairs of 2 stimuli = 6 choices
//   // values: initial value of each choice option (length nChoices)
//   // VV: output for trial-by-trial values of each choice option. In R, a matrix of size (nTrials, nChoices); here in C, an array of length (nTrials*nChoices)
//   // PE: output for trial-by-trial prediction errors for each choie option. Identical size as VV
//   // outcome: (nTrials, nChoices): trial by trial outcomes per choice option. NA is no output.
//   // eta1, eta2: floats with the learning rates (positive and negative, respectively)
//
//   // declare some variables
//   //  double updateValue;
//   double this_lr;
//   double this_variance;
//   int mat_idx = 0;
//   int nt = *nTrials;
//   int nc = *nChoices;
//   static double dv[64] = {0};   // WARNING: MAXIMUM NUMBER OF CHOICES IS FIXED HERE TO 64. Anything more will crash.
//
//   // Loop over trials i
//   for(unsigned int i = 0; i < nt; i++) {
//     //    printf("Trial N: %d\n", i);
//     //    this_eta1 = eta1[i];  // trialwise learning rate for positive PE
//     //    this_eta2 = eta2[i];  // trialwise learning rate for negative PE
//
//     // Loop over choice alternatives
//     for(unsigned int ch = 0; ch < nc; ch++) {
//       //      printf("Choice option: %d\n", ch);
//       mat_idx = nt*ch+i;  // Where in the VV-matrix and outcome-matrix are we?
//       adaptedValues[mat_idx] = values[ch]; // Offload current values
//
//       // If the outcome for this choice was not NA, update (note that o_ == o_ returns FALSE if o_ is NA)
//       if(outcome[mat_idx] == outcome[mat_idx]) {
//         this_lr = learningRates[mat_idx];
//         this_variance = values[ch] * (1-values[ch]);  // ASSUMES BINOMIAL DISTRIBUTION!
//         //        printf("Outcome is: %.3f, this is NOT NA!\n", outcome[mat_idx]);
//         dv[ch] = outcome[mat_idx] - values[ch];  // prediction error dv = outcome choice - value of choice
//         predictionErrors[mat_idx] = dv[ch];  // offload PE
//
//         // Do we update with eta1 or eta2?
//         //        updateValue = dv[ch] > 0 ? this_eta1 : this_eta2;  // Ternary expression  (if-else in a one-liner)
//
//         values[ch] = values[ch] + this_variance*this_lr*dv[ch];  // Update value
//
//         // offload actual learning rate
//         learningRates[mat_idx] = this_variance*this_lr;
//       }
//     }
//   }
// }
