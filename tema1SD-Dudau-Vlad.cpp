#include <bits/stdc++.h>

#include <chrono>

#define N 10000001

using namespace std;
using namespace std::chrono;

ifstream fin("teste.in");
ofstream fout("teste.out");
int32_t tablou[N],copie_tablou[N];
int32_t auxiliar[N];
int32_t frecventa[N];
int nr_numere, val_max;
const int RADIX = 255;
/// Bubble sort
void BubbleSort(int v[],int nr_elemente)
{
    int sortate = 0;
    while (sortate == 0)
    {
        sortate = 1;
        for (int i = 0; i < nr_elemente - 1; i++)
            if (v[i] > v[i + 1])
            {
                swap(v[i], v[i + 1]);
                sortate = 0;
            }
    }
}

/// QuickSort1

int Pivotare(int v[],int stanga, int dreapta)
{
    int i, j, pasi = 0, pasj = 1;
    i = stanga;
    j = dreapta;
    while (i < j)
    {
        if (v[i] > v[j])
        {
            swap(v[i], v[j]);
            pasi = 1 - pasi;
            pasj = 1 - pasj;
        }
        i = i + pasi;
        j = j - pasj;
    }
    return i;
}
void QuickSort1(int v[],int stanga, int dreapta)
{
    if (stanga < dreapta)
    {
        int pivot = Pivotare(v,stanga, dreapta);
        QuickSort1(v,stanga, pivot - 1);
        QuickSort1(v,pivot + 1, dreapta);
    }
}

///QuickSort cu pivot randomizat

int PartitionareRandomizata (int v[],int stanga, int dreapta )
{

    int i, j;
    i = stanga-1;
    j= dreapta+1;
    int randpoz, r1, r2, r3;
    int nr = dreapta-stanga+1;

    r1 = rand() % nr + stanga;
    r2 = rand() % nr + stanga;
    r3 = rand() % nr + stanga;

    int minr, maxr;
    minr = min ( min(r1, r2), min(r2, r3) );
    maxr = max ( max(r1, r2), max(r2, r3) );

    randpoz = r1+r2+r3 - minr-maxr;

    int pivot = v[randpoz];

    while (true)
    {
        do
        {
            ++i;
        }
        while( v[i] < pivot);
        do
        {
            --j;
        }
        while( v[j] > pivot);
        if( i>=j )
            return j;
        swap(v[i], v[j]);
    }
}
void QuickSort3(int v[],int stanga, int dreapta)
{
    if(stanga<dreapta)
    {
        int k = PartitionareRandomizata(v,stanga,dreapta);
        QuickSort3(v,stanga, k);
        QuickSort3(v,k+1, dreapta);
    }
}


/// QuickSort cu pivot mediana

int partitionare(int v[], int stanga, int dreapta, int pozitie);
int kminim(int v[], int stanga, int dreapta, int k);
void QuickSort2(int v[], int stanga, int sfarsit)
{
    if (stanga < sfarsit)
    {
        /// Lungimea subvectorului
        int lungime = sfarsit - stanga + 1;

        /// Caut mediana vectorlui
        int med = kminim(v, stanga, sfarsit, lungime / 2);

        /// Partitionez vectorul in jurul medianei
        int p = partitionare(v, stanga, sfarsit, med);

        QuickSort2(v, stanga, p - 1);
        QuickSort2(v, p + 1, sfarsit);
    }
}
/// sorteaza maxim 5 elemente de fiecare data
int CautaMediana(int v[], int n)
{
    sort(v, v + n);   /// Sorteaza vectorul
    return v[n / 2];  /// Returneaza elementul din mijloc
}

/// Returneaza al k-lea minim
int kminim(int v[], int stanga, int dreapta, int k)
{
    /// Daca k < nr elemente din vector
    if (k > 0 && k <= dreapta - stanga + 1)
    {
        int numar =
            dreapta - stanga + 1;  /// Numarul de elemente de la [l..dreapta]

        /// Impart vectorul in grupuri de cate 5 calculand mediana
        int i, mediana[(numar + 4) / 5];
        for (i = 0; i < numar / 5; i++)
            mediana[i] = CautaMediana(v + stanga + i * 5, 5);
        if (i * 5 < numar)  /// Pentru ultimul grup care poate avea mai putin de
            /// 5 elemente
        {
            mediana[i] = CautaMediana(v + stanga + i * 5, numar % 5);
            i++;
        }

        int mediana_medianelor =
            (i == 1) ? mediana[i - 1] : kminim(mediana, 0, i - 1, i / 2);

        int pos = partitionare(v, stanga, dreapta, mediana_medianelor);

        /// Daca pozitia e aceeasi cu k
        if (pos - stanga == k - 1)
            return v[pos];
        if (pos - stanga > k - 1)
            return kminim(v, stanga, pos - 1, k);

        return kminim(v, pos + 1, dreapta, k - pos + stanga - 1);
    }

    return INT_MAX;
}

void swap(int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}
int partitionare(int v[], int stanga, int dreapta, int element)
{
    int i;
    for (i = stanga; i < dreapta; i++)
        if (v[i] == element)
            break;
    swap(&v[i], &v[dreapta]);

    i = stanga;
    for (int j = stanga; j <= dreapta - 1; j++)
    {
        if (v[j] <= element)
        {
            swap(&v[i], &v[j]);
            i++;
        }
    }
    swap(&v[i], &v[dreapta]);
    return i;
}

/// MergeSort


void Interclasare(int *v, int s, int m, int d)
{
    int i, j, k;
    i = s;
    j = m + 1;
    k = 0;
    while (i <= m && j <= d)
        if (v[i] < v[j])
            auxiliar[++k] = v[i++];
        else
            auxiliar[++k] = v[j++];
    while (i <= m)
        auxiliar[++k] = v[i++];
    while (j <= d)
        auxiliar[++k] = v[j++];
    k = 0;
    for (i = s; i <= d; i++)
        v[i] = auxiliar[++k];
}
void MergeSort(int v[], int s, int d)
{
    if (s < d)
    {
        int m = (s + d) / 2;
        MergeSort(v, s, m);
        MergeSort(v, m + 1, d);
        Interclasare(v, s, m, d);
    }
}

/// CountSort

void countSort(int v[], int n)
{
    int maxi = v[0];
    for (int i = 1; i < n; i++)
    {
        if (v[i] > maxi)
            maxi = v[i];
    }

    for (int i = 0; i <= maxi; ++i)
    {
        frecventa[i] = 0;
    }

    for (int i = 0; i < n; i++)
    {
        frecventa[v[i]]++;
    }
    for (int i = 1; i <= maxi; i++)
    {
        frecventa[i] += frecventa[i - 1];
    }
    for (int i = n - 1; i >= 0; i--)
    {
        auxiliar[frecventa[v[i]] - 1] = v[i];
        frecventa[v[i]]--;
    }
    for (int i = 0; i < n; i++)
    {
        v[i] = auxiliar[i];
    }
}

/// Radix Sort avand la baza count sort

int ElementulMaxim(int v[], int n)
{
    int maxi = v[0];
    for (int i = 1; i < n; i++)
        if (v[i] > maxi)
            maxi = v[i];
    return maxi;
}

/// Counting Sort - se afla la baza RadixSort
void CountSortRadix(int v[], int n, int exp)
{
    int i, numara[32] = {0};

    /// Tinem minte numarul de aparitii in vectorul numara
    for (i = 0; i < n; i++)
        numara[(v[i] / exp) % 32]++;

    /// Modificam numara[i] astfel incat numara[i] sa contina pozitia cifrei din auxiliar[]
    for (i = 1; i < 32; i++)
        numara[i] += numara[i - 1];

    /// Creem vectorul auxiliar
    for (i = n - 1; i >= 0; i--)
    {
        auxiliar[numara[(v[i] / exp) % 32] - 1] = v[i];
        numara[(v[i] / exp) % 32]--;
    }

    /// Copiem continutul vectorului auxiliar[] in v[]
    for (i = 0; i < n; i++)
        v[i] = auxiliar[i];
}

/// Radix Sort
void radixsort(int v[], int w[], int n)
{
    /// Calculam cel mai mare element pentru a afla numarul de cifre al acestuia
    int m = ElementulMaxim(v, n);

    ///Folosim Counting Sort pentru fiecare cifra
    for (int exp = 1; m / exp > 0; exp *= 32)
        CountSortRadix(v, n, exp);
}

/// RadixSort cu operatii pe biti

void CountSortRadixBiti(int v[], int n, int byte)
{
    int i;
    int numara[256];
    int pozitie[256];
    for (i = 0; i < 256; i++)
        numara[i] = 0;
    for (int i = 0; i < n; i++)
        ++numara[(v[i] >> byte) & RADIX];

    pozitie[0] = 0;

    for (int i = 1; i < 256; i++)
        pozitie[i] = pozitie[i - 1] + numara[i - 1];

    for (int i = 0; i < n; i++)
        auxiliar[pozitie[(v[i] >> byte) & RADIX]++] = v[i];

    for (int i = 0; i < n; i++)
        v[i] = auxiliar[i];
}
void RadixSortBiti(int v[], int n)
{
    for (int i = 0; i < 32; i += 8)
        CountSortRadixBiti(v, n, i);
}

/// Functie care verifica daca un vector este sortat

bool Vector_Sortat(int v[], int n)
{
    for (int i = 0; i < n - 1; i++)
        if (v[i] > v[i + 1])
            return false;
    return true;
}
void Copiere(int v[], int w[], int n)
{
    for (int i = 0; i < n; i++)
        v[i] = w[i];
}
int getRandomNumber(int min, int max)
{
    static constexpr double fraction { 1.0 / (RAND_MAX + 1.0) };
    return min + static_cast<int>((max - min + 1) * (std::rand() * fraction));
}
int main()
{
    int teste;
    fin >> teste;
    for (int i = 1; i <= teste; i++)
    {
        fin >> nr_numere >> val_max;
        fout << nr_numere << " " << val_max << "\n";
        for (int j = 0; j < nr_numere; j++)
        {
            tablou[j] = getRandomNumber(0,val_max);
            copie_tablou[j] = tablou[j];
        }
        fout << "\n";
        auto start = high_resolution_clock::now();
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        ///BubbleSort

        if (nr_numere <= 20000)
        {
            start = high_resolution_clock::now();
            BubbleSort(tablou,nr_numere);
            stop = high_resolution_clock::now();
            duration = duration_cast<microseconds>(stop - start);
            fout << "BubbleSort  " << duration.count() / 1000 << "  "
                 << Vector_Sortat(tablou, nr_numere) << "\n";
            Copiere(tablou, copie_tablou, nr_numere);
        }
        else
        {
            fout<<"BubbleSort nu functioneaza bine pe tablouri foarte mari, ineficient."<<"\n";
        }

        ///QuickSort1

        start = high_resolution_clock::now();
        QuickSort1(tablou,0, nr_numere - 1);
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        fout << "QuickSort " << duration.count() / 1000 << " milisecunde "<< Vector_Sortat(tablou, nr_numere) << "\n";
        Copiere(tablou, copie_tablou, nr_numere);

        ///QuickSort pivot random

        start = high_resolution_clock::now();
        QuickSort3(tablou,0, nr_numere - 1);
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        fout << "QuickSort pivot randomizat " << duration.count() / 1000 << " milisecunde "<< Vector_Sortat(tablou, nr_numere) << "\n";
        Copiere(tablou, copie_tablou, nr_numere);

        ///QuickSort mediana medianelor

        if(val_max>20000)
        {
            start = high_resolution_clock::now();
            QuickSort2(tablou, 0, nr_numere - 1);
            stop = high_resolution_clock::now();
            duration = duration_cast<microseconds>(stop - start);
            fout << "Quicksort mediana-medianelor  " << duration.count() / 1000<< " milisecunde " << Vector_Sortat(tablou, nr_numere) << "\n";
            Copiere(tablou, copie_tablou, nr_numere);
        }
        else
        {
            fout<<"QuickSort mediana-medianelor functioneaza ineficient pe numere mici";
            fout<<"\n";
        }

        ///MergeSort

        start = high_resolution_clock::now();
        MergeSort(tablou, 0, nr_numere - 1);
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        fout << "MergeSort:  " << duration.count() / 1000 << " milisecunde "
             << Vector_Sortat(tablou, nr_numere) << "\n";
        Copiere(tablou, copie_tablou, nr_numere);

        ///CountSort

        if(val_max<=100000)
        {
            start = high_resolution_clock::now();
            countSort(tablou, nr_numere);
            stop = high_resolution_clock::now();
            duration = duration_cast<microseconds>(stop - start);
            fout << "CountSort  " << duration.count() / 1000 << " milisecunde  "<< Vector_Sortat(tablou, nr_numere) << "\n";
            Copiere(tablou, copie_tablou, nr_numere);
        }
        else
        {
            fout<<"CountSort functioneaza mult mai bine pe numere mici,asadar am setat limita 100000\n";
        }

        ///RadixSort

        start = high_resolution_clock::now();
        radixsort(tablou, auxiliar, nr_numere);
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        fout << "RadixSort  " << duration.count() / 1000 << " milisecunde  "
             << Vector_Sortat(tablou, nr_numere) << "\n";
        Copiere(tablou, copie_tablou, nr_numere);

        ///RadixSort biti

        start = high_resolution_clock::now();
        RadixSortBiti(tablou, nr_numere);
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        fout << "RadixSortBiti  " << duration.count() / 1000 << " milisecunde  "
             << Vector_Sortat(tablou, nr_numere) << "\n";
        Copiere(tablou, copie_tablou, nr_numere);

        ///Sort STL

        start = high_resolution_clock::now();
        sort(tablou+0,tablou+nr_numere);
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        fout << "Functia sort din STL " << duration.count() / 1000 << " milisecunde "<< Vector_Sortat(tablou, nr_numere) << "\n";
        Copiere(tablou, copie_tablou, nr_numere);
        fout << "\n";
    }
    ///Pentru vector deja sortat

    fout<<"\n\n\nDEJA SORTAT\n";
    sort(tablou+0,tablou+nr_numere);
    auto start = high_resolution_clock::now();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    ///BubbleSort

    if (nr_numere <= 20000)
    {
        start = high_resolution_clock::now();
        BubbleSort(tablou,nr_numere);
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        fout << "BubbleSort  " << duration.count() / 1000 << "  "
             << Vector_Sortat(tablou, nr_numere) << "\n";
        Copiere(tablou, copie_tablou, nr_numere);
    }
    else
    {
        fout<<"BubbleSort nu functioneaza bine pe tablouri foarte mari, ineficient."<<"\n";
    }

    ///QuickSort1

    start = high_resolution_clock::now();
    QuickSort1(tablou,0, nr_numere - 1);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    fout << "QuickSort " << duration.count() / 1000 << " milisecunde "<< Vector_Sortat(tablou, nr_numere) << "\n";
    Copiere(tablou, copie_tablou, nr_numere);

    ///QuickSort pivot random

    start = high_resolution_clock::now();
    QuickSort3(tablou,0, nr_numere - 1);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    fout << "QuickSort pivot randomizat " << duration.count() / 1000 << " milisecunde "<< Vector_Sortat(tablou, nr_numere) << "\n";
    Copiere(tablou, copie_tablou, nr_numere);

    ///QuickSort mediana medianelor

    if(val_max>20000)
    {
        start = high_resolution_clock::now();
        QuickSort2(tablou, 0, nr_numere - 1);
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        fout << "Quicksort mediana-medianelor  " << duration.count() / 1000<< " milisecunde " << Vector_Sortat(tablou, nr_numere) << "\n";
        Copiere(tablou, copie_tablou, nr_numere);
    }
    else
    {
        fout<<"QuickSort mediana-medianelor functioneaza ineficient pe numere mici";
        fout<<"\n";
    }

    ///MergeSort

    start = high_resolution_clock::now();
    MergeSort(tablou, 0, nr_numere - 1);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    fout << "MergeSort:  " << duration.count() / 1000 << " milisecunde "
         << Vector_Sortat(tablou, nr_numere) << "\n";
    Copiere(tablou, copie_tablou, nr_numere);

    ///CountSort

    if(val_max<=100000)
    {
        start = high_resolution_clock::now();
        countSort(tablou, nr_numere);
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        fout << "CountSort  " << duration.count() / 1000 << " milisecunde  "<< Vector_Sortat(tablou, nr_numere) << "\n";
        Copiere(tablou, copie_tablou, nr_numere);
    }
    else
    {
        fout<<"CountSort functioneaza mult mai bine pe numere mici,asadar am setat limita 100000\n";
    }

    ///RadixSort

    start = high_resolution_clock::now();
    radixsort(tablou, auxiliar, nr_numere);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    fout << "RadixSort  " << duration.count() / 1000 << " milisecunde  "
         << Vector_Sortat(tablou, nr_numere) << "\n";
    Copiere(tablou, copie_tablou, nr_numere);

    ///RadixSort biti

    start = high_resolution_clock::now();
    RadixSortBiti(tablou, nr_numere);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    fout << "RadixSortBiti  " << duration.count() / 1000 << " milisecunde  "
         << Vector_Sortat(tablou, nr_numere) << "\n";
    Copiere(tablou, copie_tablou, nr_numere);

    ///Sort STL

    start = high_resolution_clock::now();
    sort(tablou+0,tablou+nr_numere);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    fout << "Functia sort din STL " << duration.count() / 1000 << " milisecunde "<< Vector_Sortat(tablou, nr_numere) << "\n";
    Copiere(tablou, copie_tablou, nr_numere);
    fout << "\n";
    return 0;
}
