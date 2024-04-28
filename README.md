## Условия
- Компилятор mingw 11.0 w64

## Инструкция по запуску программы

- Из папки src собрать решение командой 
```cmd
g++ -o <название исполняемого файла> main.cpp DataParserQuick.h GeneticAlgorithm.h GeneticAlgorithm.cpp Individuals.h Population.h Population.cpp -O3 -fopenmp
```

- Выполнить команду
```cmd
<название исполняемого файла> <путь к входному файлу> <путь к выходному файлу>
```
