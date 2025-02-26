import multiprocessing
import os
import time
import argparse


def read_virus_taxids(virus_file_path):
    try:
        virus_taxids = set()
        with open(virus_file_path, 'r', encoding='utf-8') as virus_file:
            for line in virus_file:
                line = line.strip()
                if line:
                    virus_taxids.add(line)
        print(f"Прочитано {len(virus_taxids)} вирусных TaxID из {virus_file_path}")
        return virus_taxids
    except FileNotFoundError:
        print(f"Ошибка: Файл не найден по пути: {virus_file_path}")
        exit()
    except UnicodeDecodeError:
        print(f"Ошибка: Проблемы с кодировкой файла {virus_file_path}. Попробуйте указать другую кодировку.")
        exit()
    except Exception as e:
        print(f"Произошла ошибка при чтении {virus_file_path}: {e}")
        exit()


def process_chunk(file_path, virus_taxids, chunk_start, chunk_size, queue):
    results = []
    with open(file_path, 'r', encoding='utf-8') as input_file:
        input_file.seek(chunk_start)
        for i in range(chunk_size):
            line = input_file.readline()
            if not line:
                break
            line = line.strip()
            if line:
                parts = line.split('\t')
                if len(parts) == 2:
                    protein_id, taxid = parts
                    if taxid in virus_taxids:
                        results.append(f"{protein_id}\t{taxid}\n")
    queue.put(results)

def process_multithreaded(file_path, virus_taxids, num_threads):
    file_size = os.path.getsize(file_path)
    chunk_size = file_size // num_threads

    queue = multiprocessing.Queue()

    processes = []
    for i in range(num_threads):
        chunk_start = i * chunk_size
        current_chunk_size = chunk_size if i < num_threads - 1 else file_size - chunk_start
        process = multiprocessing.Process(target=process_chunk, args=(file_path, virus_taxids, chunk_start, current_chunk_size, queue))
        processes.append(process)
        process.start()

    for process in processes:
        process.join()

    all_results = []
    while not queue.empty():
        all_results.extend(queue.get())

    return all_results


def main():
    parser = argparse.ArgumentParser(description="Фильтрация protein_taxids.txt по списку вирусных TaxID.")
    parser.add_argument("input_file", help="Путь к файлу protein_taxids.txt")
    parser.add_argument("virus_file", help="Путь к файлу virus.txt")
    parser.add_argument("output_file", help="Путь к выходному файлу")
    parser.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(),
                        help="Количество рабочих потоков (по умолчанию: количество ядер CPU)")
    parser.add_argument("-b", "--benchmark", action="store_true",
                        help="Запустить в режиме бенчмарка (только измерение времени, без записи в файл)")
    parser.add_argument("-n", "--num_trials", type=int, default=3,
                        help="Количество прогонов бенчмарка (по умолчанию: 3)")

    args = parser.parse_args()

 
    virus_taxids = read_virus_taxids(args.virus_file)

    start_time = time.time()

    
    filtered_data = process_multithreaded(args.input_file, virus_taxids, args.threads)

    end_time = time.time()

    if args.benchmark:
        print(f"Бенчмарк: Время выполнения: {end_time - start_time:.2f} секунд, Потоков: {args.threads}")
    else:
        with open(args.output_file, 'w', encoding='utf-8') as filtered_file:
            filtered_file.writelines(filtered_data)
        print(f"Всего строк отфильтровано и записано в файл: {len(filtered_data)}")
        print(f"Время выполнения: {end_time - start_time:.2f} секунд")
        print(f"Результаты записаны в файл: {args.output_file}")


if __name__ == "__main__":
    main()
