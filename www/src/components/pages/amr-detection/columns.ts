import { h } from 'vue'
import type { ColumnDef } from '@tanstack/vue-table'
import { ArrowUpDown, ArrowUp, ArrowDown } from "@lucide/vue"
import { Button } from '@/components/ui/button'

export interface AmrDetectionRow {
    sample: string
    query: string
    unitLabel: string
    callType: string
    category: string
    subtype: string
    className: string
    subclass: string
    callFraction: number
    diagnosticKmers: string
    span: string
    memberCount: number
}

function sortableHeader(label: string, align?: 'right') {
    return ({ column }: { column: { toggleSorting: (asc: boolean) => void; getIsSorted: () => string | false } }) => {
        const sorted = column.getIsSorted()
        const icon = sorted === 'asc' ? ArrowUp : sorted === 'desc' ? ArrowDown : ArrowUpDown
        return h(Button, {
            variant: 'ghost',
            class: align === 'right' ? 'w-full justify-end' : undefined,
            onClick: () => column.toggleSorting(sorted === 'asc'),
        }, () => [label, h(icon, { class: 'ml-2 h-4 w-4' })])
    }
}

export const columns: ColumnDef<AmrDetectionRow>[] = [
    {
        accessorKey: 'unitLabel',
        header: sortableHeader('Report unit'),
        cell: ({ row }) => h('div', { class: 'font-medium font-mono truncate' }, row.getValue('unitLabel')),
    },
    {
        accessorKey: 'callType',
        header: sortableHeader('Call'),
    },
    {
        accessorKey: 'category',
        header: sortableHeader('Category'),
    },
    {
        accessorKey: 'subtype',
        header: sortableHeader('Subtype'),
    },
    {
        accessorKey: 'className',
        header: sortableHeader('Class'),
    },
    {
        accessorKey: 'subclass',
        header: sortableHeader('Subclass'),
    },
    {
        accessorKey: 'callFraction',
        header: sortableHeader('Fraction', 'right'),
        cell: ({ row }) => {
            const fraction = row.getValue('callFraction') as number
            return h('div', { class: 'text-right font-medium' }, `${(fraction * 100).toFixed(1)}%`)
        },
    },
    {
        accessorKey: 'diagnosticKmers',
        header: sortableHeader('Diagnostic k-mers', 'right'),
        cell: ({ row }) => h('div', { class: 'text-right font-mono' }, row.getValue('diagnosticKmers')),
    },
    {
        accessorKey: 'query',
        header: sortableHeader('Query'),
        cell: ({ row }) => h('div', { class: 'font-mono truncate' }, row.getValue('query')),
    },
    {
        accessorKey: 'span',
        header: sortableHeader('Span', 'right'),
        cell: ({ row }) => h('div', { class: 'text-right font-mono' }, row.getValue('span')),
    },
    {
        accessorKey: 'memberCount',
        header: sortableHeader('Members', 'right'),
        cell: ({ row }) => h('div', { class: 'text-right' }, row.getValue('memberCount')),
    },
]
