package duckdb

import (
	"os"
	"time"
)

// FileFingerprint holds stat-based identity for a file.
type FileFingerprint struct {
	Path    string
	Size    int64
	ModTime time.Time
}

// StatFile creates a FileFingerprint from an on-disk file.
func StatFile(path string) (FileFingerprint, error) {
	info, err := os.Stat(path)
	if err != nil {
		return FileFingerprint{}, err
	}
	return FileFingerprint{
		Path:    path,
		Size:    info.Size(),
		ModTime: info.ModTime(),
	}, nil
}
