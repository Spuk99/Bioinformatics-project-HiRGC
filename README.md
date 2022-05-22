# Bioinformatics 1

This is a project from the 'Bioinformatics 1' course at the Faculty of Electrical Engineering and Computing ([FER](https://www.fer.unizg.hr/en)).

Course information can be found [here](https://www.fer.unizg.hr/en/course/bio1).

Goal of this project is to implement the high-speed and high-ratio referential genome compression algorithm. Original article can be found [here](https://doi.org/10.1093/bioinformatics/btx412) and the original repository can be found [here](https://github.com/yuansliu/HiRGC).

# Authors

- [Katarina Mišura](https://github.com/Spuk99)
- [Marko Marfat](https://github.com/mmarfat)

# Necessary prerequisites 

- [MinGW](https://www.mingw-w64.org/downloads/)
- [7zip](https://www.7-zip.org/)

<pre>
<b>NOTE:</b> This code is meant to be run on a Linux operating system.
</pre>

# Compile

```
g++ hirgc.cpp -o hirgc -O3
g++ de_hirgc.cpp -o de_hirgc -O3
```

# Run

<pre>
<b>NOTE:</b> Use the provided test2ref.fa as reference and test2tar.fa as target.
</pre>

```
./hirgc -r reference.fa -t target.fa
./de_hirgc -r reference.fa -t target.7z
```


