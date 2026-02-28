package main

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/spf13/cobra"
	"github.com/spf13/viper"
	"gopkg.in/yaml.v3"
)

func newConfigCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "config",
		Short: "Manage vibe-vep configuration",
		Long:  "Show, get, or set configuration values. Config is stored in ~/.vibe-vep.yaml.",
		Example: `  vibe-vep config                            # show all config
  vibe-vep config set annotations.alphamissense true  # enable AlphaMissense
  vibe-vep config get annotations.alphamissense       # get a value`,
		Args: cobra.NoArgs,
		RunE: func(cmd *cobra.Command, args []string) error {
			return runConfigShow()
		},
	}

	cmd.AddCommand(newConfigSetCmd())
	cmd.AddCommand(newConfigGetCmd())

	return cmd
}

func newConfigSetCmd() *cobra.Command {
	return &cobra.Command{
		Use:   "set <key> <value>",
		Short: "Set a configuration value",
		Args:  cobra.ExactArgs(2),
		RunE: func(cmd *cobra.Command, args []string) error {
			return runConfigSet(args[0], args[1])
		},
	}
}

func newConfigGetCmd() *cobra.Command {
	return &cobra.Command{
		Use:   "get <key>",
		Short: "Get a configuration value",
		Args:  cobra.ExactArgs(1),
		RunE: func(cmd *cobra.Command, args []string) error {
			return runConfigGet(args[0])
		},
	}
}

func runConfigShow() error {
	settings := viper.AllSettings()
	if len(settings) == 0 {
		fmt.Println("# No configuration set. Config file: ~/.vibe-vep.yaml")
		return nil
	}

	out, err := yaml.Marshal(settings)
	if err != nil {
		return fmt.Errorf("marshaling config: %w", err)
	}
	fmt.Print(string(out))
	return nil
}

func runConfigSet(key, value string) error {
	// Parse boolean-like values
	switch value {
	case "true", "yes", "on":
		viper.Set(key, true)
	case "false", "no", "off":
		viper.Set(key, false)
	default:
		viper.Set(key, value)
	}

	// Ensure config file exists
	cfgFile := viper.ConfigFileUsed()
	if cfgFile == "" {
		home, err := os.UserHomeDir()
		if err != nil {
			return fmt.Errorf("cannot determine home directory: %w", err)
		}
		cfgFile = filepath.Join(home, ".vibe-vep.yaml")
	}

	if err := viper.WriteConfigAs(cfgFile); err != nil {
		return fmt.Errorf("writing config: %w", err)
	}

	fmt.Printf("Set %s = %s in %s\n", key, value, cfgFile)
	return nil
}

func runConfigGet(key string) error {
	val := viper.Get(key)
	if val == nil {
		return fmt.Errorf("key %q is not set", key)
	}
	fmt.Println(val)
	return nil
}
