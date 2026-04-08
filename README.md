## Структура проекта

Проект разделён на модули: геометрия и сетка, акустика, динамика пузырька, связка между ними, ввод-вывод и постобработка.

### Верхний уровень

- `include/` — заголовочные файлы (`.hpp`) с интерфейсами классов, структур и функций.
- `src/` — реализации (`.cpp`) для файлов из `include/`.
- `apps/` — отдельные исполняемые программы для запуска частей проекта.

---

### `include/core/` и `src/core/`

Базовые служебные части проекта.

- `types.hpp` — общие типы и алиасы, используемые по всему проекту.
- `constants.hpp` — физические и численные константы по умолчанию.
- `config.hpp` / `config.cpp` — структуры конфигурации и работа с параметрами запуска.
- `logger.hpp` / `logger.cpp` — вывод логов и служебных сообщений.

---

### `include/geometry/` и `src/geometry/`

Описание геометрии сосуда и отражателя.

- `geometry_builder.hpp` / `geometry_builder.cpp` — общий сборщик геометрии задачи.
- `vessel_geometry.hpp` / `vessel_geometry.cpp` — описание формы сосуда.
- `reflector_geometry.hpp` / `reflector_geometry.cpp` — описание параболического отражателя.

---

### `include/mesh/` и `src/mesh/`

Работа с сеткой и разметкой границ.

- `gmsh_driver.hpp` / `gmsh_driver.cpp` — генерация сетки через Gmsh.
- `mesh_tags.hpp` / `mesh_tags.cpp` — физические метки границ и областей (стенка, источник, отражатель и т.д.).

---

### `include/acoustics/` и `src/acoustics/`

Акустическая задача в жидкости.

- `acoustics_problem.hpp` / `acoustics_problem.cpp` — полная постановка акустической задачи.
- `helmholtz_solver.hpp` / `helmholtz_solver.cpp` — решение уравнения Гельмгольца.
- `boundary_conditions.hpp` / `boundary_conditions.cpp` — граничные условия на стенках, источнике и отражателе.
- `field_sampler.hpp` / `field_sampler.cpp` — извлечение значений поля давления в заданных точках.
- `resonance_scan.hpp` / `resonance_scan.cpp` — сканирование по частоте и поиск резонансных режимов.

---

### `include/bubble/` и `src/bubble/`

Модель динамики пузырька.

- `rayleigh_plesset.hpp` / `rayleigh_plesset.cpp` — уравнение Рэя–Плессета.
- `gas_model.hpp` / `gas_model.cpp` — модель газа внутри пузырька.
- `thermal_model.hpp` / `thermal_model.cpp` — теплообмен и температурная модель газа.
- `bubble_rhs.hpp` / `bubble_rhs.cpp` — правая часть системы ОДУ для пузырька.
- `ode_solver.hpp` / `ode_solver.cpp` — численное решение системы ОДУ.

---

### `include/coupling/` и `src/coupling/`

Связка акустики и пузырька.

- `pressure_coupling.hpp` / `pressure_coupling.cpp` — передача акустического давления в модель пузырька.
- `pipeline.hpp` / `pipeline.cpp` — полный вычислительный сценарий: геометрия → акустика → пузырёк.

---

### `include/io/` и `src/io/`

Ввод-вывод данных.

- `yaml_reader.hpp` / `yaml_reader.cpp` — чтение параметров из YAML-конфига.
- `csv_writer.hpp` / `csv_writer.cpp` — сохранение численных результатов в CSV.
- `vtk_writer.hpp` / `vtk_writer.cpp` — экспорт полей и данных для визуализации.
- `result_writer.hpp` / `result_writer.cpp` — сохранение итоговых результатов запуска.

---

### `apps/`

Отдельные точки входа в проект.

- `build_geometry.cpp` — только построение геометрии и сетки.
- `solve_acoustics.cpp` — только решение акустической задачи.
- `solve_bubble.cpp` — только решение модели пузырька.
- `run_pipeline.cpp` — полный запуск всей связанной модели.
- `scan_frequency.cpp` — сканирование по частоте.

### Логика проекта

Проект устроен по цепочке:

1. В `geometry/` и `mesh/` строится область расчёта.
2. В `acoustics/` решается акустическая задача в жидкости.
3. В `coupling/` давление из акустики передаётся в модель пузырька.
4. В `bubble/` решается система ОДУ для динамики пузырька и газа.
5. В `io/` сохраняются результаты.
