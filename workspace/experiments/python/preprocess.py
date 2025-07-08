import argparse
import json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot configuration from JSON file")
    parser.add_argument(
        "--base_config_file", type=str, help="Path to the base configuration JSON file"
    )
    parser.add_argument(
        "--task_config_file", type=str, help="Path to the task configuration JSON file"
    )
    args = parser.parse_args()

    # Load the configuration file
    with open(args.base_config_file, "r") as f:
        base_config = json.load(f)
    with open(args.task_config_file, "r") as f:
        config = json.load(f)

    # Merge the base configuration with the task configuration
    config = {**base_config, "tasks": config.get("tasks", {})}
    save_path = args.task_config_file
    with open(save_path, "w") as f:
        json.dump(config, f, indent=4)
