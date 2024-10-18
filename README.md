Code for fitting and analyzing multi-scale racing diffusion models (MS-RDMs)

Compared to the paper, naming conventions here are different:
- Accuracy memory (AM) in the paper is named accuracy history (AH) in the code
- Fluency memory (FM) in the paper is named RS in the code
- db (threshold difference) in the paper is named z here (c.f. starting point in a DDM)
- dr (rate difference) in the paper is named v here (c.f. drift rate in DDM)
- b (average threshold) in the paper is named a here (c.f. the threshold in DDM)

Abbrevations in the code roughly combine the EAM parameter that is affected, with the mechanism. I.e.,
- zSM means that SM affects the z parameter (threshold bias)
- vSM means that SM affects the v parameter (rate bias)
- uAH means that AM affects the urgency parameter u
- aAH means that AM affects the threshold parameter a
- vAH means that AM affects the rate quality parameter
- bV means that FM affects overall threshold b (apologies for the naming inconsistencies!)

Trends/cosine models are interchangeably called 'DCT' models. They're NULL if no cosines/trends were estimated.
