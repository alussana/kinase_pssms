## Kinase domain HMM

* 452/634 kinases matched in at least one subsequence

  * 38 kinases have multiple matched kinase domain (number of matches, AC)

    ```
          3 O00311
          2 O60674
          2 O75582
          2 O75676
          2 O95835
          2 P07333
          2 P09619
          2 P10721
          2 P16234
          2 P17948
          2 P23458
          2 P29597
          2 P35626
          2 P35916
          2 P35968
          2 P36888
          2 P51812
          2 P52333
          2 P78362
          2 Q13237
          2 Q15349
          2 Q15418
          2 Q15772
          2 Q56UN5
          2 Q5VST9
          2 Q8TF76
          2 Q96GX5
          2 Q96S38
          2 Q96SB4
          2 Q9BQI3
          2 Q9BUB5
          2 Q9H5K3
          2 Q9NRM7
          2 Q9NZJ5
          3 Q9P2K8
          2 Q9UK32
          2 Q9UPE1
          2 Q9Y6S9
    ```

    

* Heuristic A-loop search results in 342/452 matches, covering 335 unique kinases

  * the 7 kinases that have multiple a-loops in multiple kinase domains are:

    ```
    >O75676_33-301
    DFGLSKEFLTEEKERTFSFCGTIEYMAPE
    >O75676_416-673
    DFGFARLRPQSPGVPMQTPCFTLQYAAPE
    --
    >P51812_68-327
    DFGLSKESIDHEKKAYSFCGTVEYMAPE
    >P51812_422-679
    DFGFAKQLRAENGLLMTPCYTANFVAPE
    --
    >Q15349_59-318
    DFGLSKEAIDHDKRAYSFCGTIEYMAPE
    >Q15349_415-672
    DFGFAKQLRAGNGLLMTPCYTANFVAPE
    --
    >Q15418_62-320
    DFGLSKEAIDHEKKAYSFCGTVEYMAPE
    >Q15418_418-675
    DFGFAKQLRAENGLLMTPCYTANFVAPE
    --
    >Q15772_1601-1854
    DFGNAQELTPGEPQYCQYGTPEFVAPE
    >Q15772_2974-3218
    DFGSAQPYNPQALRPLGHRTGTLEFMAPE
    --
    >Q5VST9_6468-6721
    DFGFAQNITPAELQFSQYGSPE
    >Q5VST9_7675-7924
    DLGNAQSLSQEKVLPSDKFKDYLETMAPE
    --
    >Q9UK32_73-330
    DFGLSKESVDQEKKAYSFCGTVEYMAPE
    >Q9UK32_426-683
    DFGFAKQLRGENGLLLTPCYTANFVAPE
    ```

  * those kinases' a-loops won't be considered due to ambiguity

## UniProt annotation

* 444/634 kinases with annotation found
* Heuristic A-loop search results in 324/444 matches and unique kinases

## Overlap between HMMER-based and UniProt A-loops

UniProt ACs where HMMER-based A-loop has been retrieved but not UniProt-based. All of these A-loops are nice and short. UniProt didn't have annotation for a kinase domain for these

```
 Entry                                           Sequence Kinase domain A-loop                                      A-loop HMM
O60674  MGMACLTMTEMEGTSTSSIYQNGDISGNANSMKQIDPVLQVYLYHS...           NaN    NaN                 DFGLTKVLPQDKEYYKVKEPGESPIFWYAPE
O75582  MEEEGGSSGGAAGTSADGGDGGEQLLTVKHELRTANLTGHAEKVGI...           NaN    NaN                    DFGFARLKPPDNQPLKTPCFTLHYAAPE
P23458  MQYLNIKEDCNAMAFCAKMRSSKKTEVNLEAPEPGVEVIFYLSDRE...           NaN    NaN                 DFGLTKAIETDKEYYTVKDDRDSPVFWYAPE
P29597  MPLRHWGMARGSKPVGDGAQPMAAMGGLKVLLHWAGPGGGEPWVTF...           NaN    NaN                 DFGLAKAVPEGHEYYRVREDGDSPVFWYAPE
P52333  MAPPSEETPLIPQRSCSLLSTEAGALHVLLPARGPGPPQRLSFSFG...           NaN    NaN                 DFGLAKLLPLDKDYYVVREPGQSPIFWYAPE
Q9P2K8  MAGGRGAPGRGRDEPPESYPQRQDHELQALEAIYGADFQDLRPDAC...           NaN    NaN  DFGLATDHLAFSADSKQDDQTGDLIKSDPSGHLTGMVGTALYVSPE
```

UniProt ACs where where HMMER-based A-loop has not been retrieved, while in UniProt annotation they could be extracted. These are only two kinases (LATS1 and LATS2), and the putative A-loop is 75 AA long

```
 Entry                                           Sequence                                      Kinase domain                                             A-loop A-loop HMM
O95835  MKRSEKPEGYRQMRPKTFPASNYTVSSRQMLQEIRESLRNLSKPSD...  FVKIKTLGIGAFGEVCLARKVDTKALYATKTLRKKDVLLRNQVAHV...  DFGLCTGFRWTHDSKYYQSGDHPRQDSMDFSNEWGDPSSCRCGDRL...        NaN
Q9NRM7  MRPKTFPATTYSGNSRQRLQEIREGLKQPSKSSVQGLPAGPNSDTS...  FVKIKTLGIGAFGEVCLACKVDTHALYAMKTLRKKDVLNRNQVAHV...  DFGLCTGFRWTHNSKYYQKGSHVRQDSMEPSDLWDDVSNCRCGDRL...        NaN
```

