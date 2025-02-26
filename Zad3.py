import os
import time

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

def process_singlethreaded(file_path, virus_taxids):
    results = []
    try:
        with open(file_path, 'r', encoding='utf-8') as input_file:
            for line in input_file:
                line = line.strip()
                if line:
                    parts = line.split('\t')
                    if len(parts) == 2:
                        protein_id, taxid = parts
                        if taxid in virus_taxids:
                            results.append(f"{protein_id}\t{taxid}\n")
    except FileNotFoundError:
        print(f"Ошибка: Файл не найден по пути: {file_path}")
        exit()
    except UnicodeDecodeError:
        print(f"Ошибка: Проблемы с кодировкой файла {file_path}. Попробуйте указать другую кодировку.")
        exit()
    except Exception as e:
        print(f"Произошла ошибка при чтении {file_path}: {e}")
        exit()

    return results


def main():
    input_file_path = "protein_ids.txt"
    virus_file_path = "virus.txt"
    output_file_path = "filtered_protein_taxids.txt"

    virus_taxids = read_virus_taxids(virus_file_path)

    start_time = time.time()

    filtered_data = process_singlethreaded(input_file_path, virus_taxids)

    end_time = time.time()

    try:
        with open(output_file_path, 'w', encoding='utf-8') as filtered_file:
            filtered_file.writelines(filtered_data)
        print(f"Всего строк отфильтровано и записано в файл: {len(filtered_data)}")
        print(f"Время выполнения: {end_time - start_time:.2f} секунд")
        print(f"Результаты записаны в файл: {output_file_path}")
    except Exception as e:
        print(f"Ошибка при записи в файл: {e}")
        exit()


if __name__ == "__main__":
    main()
